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

#ifndef OB_IO_BINARY_INCLUDED
#define OB_IO_BINARY_INCLUDED

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <vector>
#include <string>

#ifndef __sgi
#include <iostream>
#else
#include <iostream.h>
#endif

namespace OpenBabel {

/* generic binary readers for basic data types */
unsigned int OB_io_read_binary(char *ccc, char *x, unsigned int size, unsigned int count);
unsigned int OB_io_read_binary(std::istream &, char *x, unsigned int size, unsigned int count);

/* generic binary writers for basic data types */
unsigned int OB_io_write_binary(char *ccc, const char *x, unsigned int size, unsigned int count);
unsigned int OB_io_write_binary(std::ostream &, const char *x, unsigned int size, unsigned int count);

/* binary readers for STL strings */
unsigned int OB_io_read_binary(char* ccc, std::string& str);
unsigned int OB_io_read_binary(std::istream& ifs, std::string& str);

/* binary writers for STL strings */
unsigned int OB_io_write_binary(char* ccc, const std::string& str);
unsigned int OB_io_write_binary(std::ostream& ofs, const std::string& str);

/* binary writers for compressed data */
unsigned int OB_io_write_binary_compressed(char*ccc, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Writes 1+(bits*Nvalues)/8 bytes
unsigned int OB_io_write_binary_compressed(std::ostream&, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Writes 1+(bits*Nvalues)/8 bytes
unsigned int OB_io_write_binary_compressed(char*ccc, float *x, unsigned int bits, unsigned int Nvalues); //Writes 9+(bits*Nvalues)/8 bytes
unsigned int OB_io_write_binary_compressed(std::ostream&, float *x, unsigned int bits, unsigned int Nvalues); //Writes 9+(bits*Nvalues)/8 bytes

/* binary readers for compressed data */
unsigned int OB_io_read_binary_compressed(char*ccc, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Reads 1+(bits*Nvalues)/8 bytes
unsigned int OB_io_read_binary_compressed(std::istream&, unsigned int *x, unsigned int bits, unsigned int Nvalues); //Reads 1+(bits*Nvalues)/8 bytes
unsigned int OB_io_read_binary_compressed(char*ccc, float *x, unsigned int bits, unsigned int Nvalues); //Reads 9+(bits*Nvalues)/8 bytes
unsigned int OB_io_read_binary_compressed(std::istream&, float *x, unsigned int bits, unsigned int Nvalues); //Reads 9+(bits*Nvalues)/8 bytes

/* utilities for reading/writing compressed data */
unsigned int OB_io_util_calc_NumBits(unsigned int *x, unsigned int N);
unsigned int OB_io_util_calc_NumBits(float *x, unsigned int N, float res);

/*!
**\brief An abstract base class with binary read and write function.
**Members of this class understand how to read and write themselves
**to and from platform independent binary streams and arrays using
**the functions of this class.
*/
class OBBinaryIO {
    public:
        virtual ~OBBinaryIO() {}
        virtual unsigned int WriteBinary(char *ccc) const = 0;
        virtual unsigned int ReadBinary(char *ccc) = 0;
        virtual unsigned int WriteBinary(std::ostream& ostr) const = 0;
        virtual unsigned int ReadBinary(std::istream& istr) = 0;
        virtual unsigned int BinarySize() const = 0;
  };

/*!
**\brief Writes a OBBinaryIO object to a file stream, preappending
**the record with the size of the OBBinaryIO objects binary record.
*/
std::ostream& operator<<(std::ostream& ostr, const OBBinaryIO& obj);

/*!
**\brief Reads in a OBBinaryIO object that has the size of the 
**binary record preappended to the record.
*/
std::istream& operator>>(std::istream& istr, OBBinaryIO& obj);

}//End namespace OpenBabel
#endif

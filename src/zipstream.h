/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
*/

/* Modified 2005 by Geoffrey R. Hutchison for modern compiliation in GCC */

#ifndef ZIPSTREAM_H
#define ZIPSTREAM_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <zlib.h>

namespace zlib_stream{

/// default gzip buffer size,
/// change this to suite your needs
const size_t default_buffer_size = 4096;

/// Compression strategy, see zlib doc.
enum EStrategy
{
	StrategyFiltered = 1,
	StrategyHuffmanOnly = 2,
	DefaultStrategy = 0
};

/** \brief A stream decorator that takes raw input and zips it to a ostream.

The class wraps up the inflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_streambuf : public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
    typedef Elem char_type;
	typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
	typedef std::vector<char_type, char_allocator_type > char_vector_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
    basic_zip_streambuf(
		ostream_reference ostream_,
		size_t level_,
		EStrategy strategy_,
		size_t window_size_,
		size_t memory_level_,
		size_t buffer_size_
		);
	
	~basic_zip_streambuf();

	int sync ();
    int overflow (int c);

	/** flushes the zip buffer and output buffer.

	This method should be called at the end of the compression. Calling flush multiple times, will lower the
	compression ratio.
	*/
	std::streamsize flush();
	/// returns a reference to the output stream
	ostream_reference get_ostream() const	{	return m_ostream;};
	/// returns the latest zlib error status
	int get_zerr() const					{	return m_err;};
	/// returns the crc of the input data compressed so far.
	long get_crc() const					{	return m_crc;};
	/// returns the size (bytes) of the input data compressed so far.
	long get_in_size() const				{	return m_zip_stream.total_in;};
	/// returns the size (bytes) of the compressed data so far.
	long get_out_size() const				{	return m_zip_stream.total_out;};
private:
	bool zip_to_stream( char_type*, std::streamsize);
	size_t fill_input_buffer();

	ostream_reference m_ostream;
	z_stream m_zip_stream;
    int m_err;
	byte_vector_type m_output_buffer;
	char_vector_type m_buffer; 
	long m_crc;
};

/** \brief A stream decorator that takes compressed input and unzips it to a istream.

The class wraps up the deflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_unzip_streambuf : 
	public std::basic_streambuf<Elem, Tr> 
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
	typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
	typedef byte_type* byte_buffer_type;
    typedef Elem char_type;
	typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
	typedef std::vector<char_type, char_allocator_type > char_vector_type;

     /** Construct a unzip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
	 basic_unzip_streambuf(
		istream_reference istream_,
		size_t window_size_,
		size_t read_buffer_size_,
		size_t input_buffer_size_
		);
	
	~basic_unzip_streambuf();

    int underflow();


	/// returns the compressed input istream
	istream_reference get_istream()	{	return m_istream;};
	/// returns the zlib stream structure
	z_stream& get_zip_stream()		{	return m_zip_stream;};
	/// returns the latest zlib error state
	int get_zerr() const					{	return m_err;};
	/// returns the crc of the uncompressed data so far 
	long get_crc() const					{	return m_crc;};
	/// returns the number of uncompressed bytes
	long get_out_size() const				{	return m_zip_stream.total_out;};
	/// returns the number of read compressed bytes
	long get_in_size() const				{	return m_zip_stream.total_in;};
private:
	void put_back_from_zip_stream();
	std::streamsize unzip_from_stream( char_type*, std::streamsize);

	size_t fill_input_buffer();

	istream_reference m_istream;
	z_stream m_zip_stream;
    int m_err;
	byte_vector_type m_input_buffer;
	char_vector_type m_buffer; 
	long m_crc;
};

/*! \brief Base class for zip ostreams

Contains a basic_zip_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
	typedef basic_zip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > zip_streambuf_type;

    /** Construct a zip stream
     * More info on the following parameters can be found in the zlib documentation.
     */
	basic_zip_ostreambase( 
		ostream_reference ostream_,
		size_t level_,
		EStrategy strategy_,
		size_t window_size_,
		size_t memory_level_,
		size_t buffer_size_
		)
		: m_buf(ostream_,level_,strategy_,window_size_,memory_level_,buffer_size_)
	{
		init(&m_buf );
	};
	
	/// returns the underlying zip ostream object
	zip_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the zlib error state
	int get_zerr() const					{	return m_buf.get_err();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the compressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the uncompressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	zip_streambuf_type m_buf;
};

/*! \brief Base class for unzip istreams

Contains a basic_unzip_streambuf.
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
	typedef std::basic_istream<Elem, Tr>& istream_reference;
	typedef basic_unzip_streambuf<
        Elem,
        Tr,
        ElemA,
        ByteT,
        ByteAT
        > unzip_streambuf_type;

	basic_zip_istreambase( 
		istream_reference istream_,
		size_t window_size_,
		size_t read_buffer_size_,
		size_t input_buffer_size_
		)
		: m_buf(istream_,window_size_, read_buffer_size_, input_buffer_size_)
	{
		init(&m_buf );
	};
	
	/// returns the underlying unzip istream object
	unzip_streambuf_type* rdbuf() { return &m_buf; };

	/// returns the zlib error state
	int get_zerr() const					{	return m_buf.get_zerr();};
	/// returns the uncompressed data crc
	long get_crc() const					{	return m_buf.get_crc();};
	/// returns the uncompressed data size
	long get_out_size() const				{	return m_buf.get_out_size();};
	/// returns the compressed data size
	long get_in_size() const				{	return m_buf.get_in_size();};
private:
	unzip_streambuf_type m_buf;
};

/*! \brief A zipper ostream

This class is a ostream decorator that behaves 'almost' like any other ostream.

At construction, it takes any ostream that shall be used to output of the compressed data.

When finished, you need to call the special method zflush or call the destructor 
to flush all the intermidiate streams.

Example:
\code
// creating the target zip string, could be a fstream
ostringstream ostringstream_;
// creating the zip layer
zip_ostream zipper(ostringstream_);

	
// writing data	
zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// zip ostream needs special flushing...
zipper.zflush();
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>	
class basic_zip_ostream : 
	public basic_zip_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_ostream<Elem,Tr>
{
public:
	typedef basic_zip_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> zip_ostreambase_type;
	typedef std::basic_ostream<Elem,Tr> ostream_type;
	typedef std::basic_ostream<Elem, Tr>& ostream_reference;
        typedef Elem char_type;

	/** Constructs a zipper ostream decorator
	 *
	 * \param ostream_ ostream where the compressed output is written
	 * \param is_gzip_ true if gzip header and footer have to be added
	 * \param level_ level of compression 0, bad and fast, 9, good and slower,
	 * \param strategy_ compression strategy
	 * \param window_size_ see zlib doc
	 * \param memory_level_ see zlib doc
	 * \param buffer_size_ the buffer size used to zip data

	 When is_gzip_ is true, a gzip header and footer is automatically added.
	 */
	basic_zip_ostream(
		ostream_reference ostream_, 
		int open_mode = std::ios::out, 
		bool is_gzip_ = false,
		size_t level_ = Z_DEFAULT_COMPRESSION,
		EStrategy strategy_ = DefaultStrategy,
		size_t window_size_ = 15,
		size_t memory_level_ = 8,
		size_t buffer_size_ = default_buffer_size
		)
	: 
		zip_ostreambase_type(
            ostream_,
            level_,
            strategy_,
            window_size_,
            memory_level_,
            buffer_size_
            ), 
		m_is_gzip(is_gzip_),
		ostream_type(this->rdbuf())
	{
		if (m_is_gzip)
			add_header();
	};
	~basic_zip_ostream()
	{
		if (m_is_gzip)
			add_footer();
	}

	/// returns true if it is a gzip 
	bool is_gzip() const		{	return m_is_gzip;};
	/// flush inner buffer and zipper buffer
	basic_zip_ostream<Elem,Tr>& zflush()	
	{	
		this->flush(); this->rdbuf()->flush(); return *this; 
	};

private:
    static void put_long(ostream_reference out_, unsigned long x_);

	void add_header();
	void add_footer();
	bool m_is_gzip;
};

/*! \brief A zipper istream

This class is a istream decorator that behaves 'almost' like any other ostream.

At construction, it takes any istream that shall be used to input of the compressed data.

Simlpe example:
\code
// create a stream on zip string
istringstream istringstream_( ostringstream_.str());
// create unzipper istream
zip_istream unzipper( istringstream_);

// read and unzip
unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;
\endcode
*/
template<
	typename Elem, 
	typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = unsigned char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_zip_istream : 
	public basic_zip_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>, 
	public std::basic_istream<Elem,Tr>
{
public:
	typedef basic_zip_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> zip_istreambase_type;
	typedef std::basic_istream<Elem,Tr> istream_type;
	typedef unsigned char byte_type;
	typedef std::basic_istream<Elem, Tr>& istream_reference;
        typedef Elem char_type;

	/** Construct a unzipper stream
	 *
	 * \param istream_ input buffer
	 * \param window_size_ 
	 * \param read_buffer_size_ 
	 * \param input_buffer_size_ 
	 */
	basic_zip_istream(
		istream_reference istream_, 
		size_t window_size_ = 15,
		size_t read_buffer_size_ = default_buffer_size,
		size_t input_buffer_size_ = default_buffer_size
		)
	  : 
		zip_istreambase_type(istream_,window_size_, read_buffer_size_, input_buffer_size_), 
		istream_type(this->rdbuf()),
		m_is_gzip(false),
		m_gzip_crc(0),
		m_gzip_data_size(0)
	{
 	      if (this->rdbuf()->get_zerr()==Z_OK)
			  check_header();
	};

	/// returns true if it is a gzip file
	bool is_gzip() const				{	return m_is_gzip;};
	/// reads the gzip header
	void read_footer();
	/** return crc check result

	When you have finished reading the compressed data, call read_footer to read the uncompressed data crc.
	This method compares it to the crc of the uncompressed data.

	\return true if crc check is succesful 
	*/
	bool check_crc() const				{	return this->get_crc() == m_gzip_crc;};
	/// return data size check
	bool check_data_size() const		{	return this->get_out_size() == m_gzip_data_size;};

	/// return the crc value in the file
	long get_gzip_crc() const			{	return m_gzip_crc;};
	/// return the data size in the file 
	long get_gzip_data_size() const		{	return m_gzip_data_size;};
protected:
    static void read_long(istream_reference in_, unsigned long& x_);

	int check_header();
	bool m_is_gzip;
	unsigned long m_gzip_crc;
	unsigned long m_gzip_data_size;
};

/// A typedef for basic_zip_ostream<char>
typedef basic_zip_ostream<char> zip_ostream;
/// A typedef for basic_zip_ostream<wchar_t>
typedef basic_zip_ostream<wchar_t> zip_wostream;
/// A typedef for basic_zip_istream<char>
typedef basic_zip_istream<char> zip_istream;
/// A typedef for basic_zip_istream<wchar_t>
typedef basic_zip_istream<wchar_t> zip_wistream;

}; // zlib_stream

#if defined(MSDOS) || (defined(WINDOWS) && !defined(WIN32))
#  define OS_CODE  0x00
#endif

#ifdef AMIGA
#  define OS_CODE  0x01
#endif

#if defined(VAXC) || defined(VMS)
#  define OS_CODE  0x02
#endif

#if defined(ATARI) || defined(atarist)
#  define OS_CODE  0x05
#endif

#ifdef OS2
#  define OS_CODE  0x06
#endif

#if defined(MACOS) || defined(TARGET_OS_MAC)
#  define OS_CODE  0x07
#endif

#ifdef TOPS20
#  define OS_CODE  0x0a
#endif

#ifdef WIN32
#  ifndef __CYGWIN__  /* Cygwin is Unix, not Win32 */
#    define OS_CODE  0x0b
#  endif
#endif

#ifdef __50SERIES /* Prime/PRIMOS */
#  define OS_CODE  0x0f
#endif

#ifndef OS_CODE
#  define OS_CODE  0x03  /* assume Unix */
#endif

#endif // HEADER


/**********************************************************************
lineend.h - Stream buffer for filtering line endings, converting \r or \r\n -> \n

 Copyright (C) 1998 by James Kanze
 Copyright (C) 2007 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_LINEEND_H
#define OB_LINEEND_H

#include <streambuf>
#include <climits>

#ifndef OBCONV
#define OBCONV
#endif

namespace OpenBabel
{

/*! \class FilteringInputStreambuf lineend.h <openbabel/lineend.h>
  \brief Delivers characters from an istream or streambuf from a source
  while filtering

  Based on an article by James Kanze, "Filtering Streambufs"
  http://kanze.james.neuf.fr/articles/fltrsbf1.html

  A FilteringInputStreambuf delivers characters on request to an istream
  or a destination rdbuf(). It receives them from a source rdbuf.
  In doing the transfer it filters them in a way decided by the class
  specified in template parameter Extractor.

  seekg and tellg requests from the stream are passed through to source rdbuf.
  This allows return to a position in the input data that was previously noted.
  This is adequate to allow OpenBabel's fastsearch indexing, but may
  not be good enough for some other applications that use random access.

  A class LineEndingExtractor converts DOS and MAC line endings to the
  UNIX line ending.

  This filtering process is potentially extendable, with a chain of
  FilteringInputStreambufs each carrying out its filtering task.
  For instance a decompression streambuf could feed a LineEnding filter,
  which in tern was read by an input stream.
*/
  template< class Extractor >
  class FilteringInputStreambuf : public std::streambuf
  {
  public:
    FilteringInputStreambuf(
      std::istream*        source = nullptr) ;
    virtual                 ~FilteringInputStreambuf()
    {
      //sync(); comment out so can be deleted in OBConversion destructor
    };
    virtual int              overflow( int ) {return EOF;};
    virtual int              underflow() ;
    virtual int              sync() ;

    //Pass the random acess functions to the source rdbuf and synchronize
    virtual std::streampos   seekoff(std::streamoff off, std::ios_base::seekdir way,
      std::ios_base::openmode which = std::ios_base::in | std::ios_base::out )
    {
      setg(  &myBuffer , &myBuffer , &myBuffer  ) ; //ensure next character is from new position
      mySource->seekg(off, way);
      std::streampos ret = mySource->tellg();
//      sync();
      return ret;
    };

    virtual std::streampos   seekpos(std::streampos sp,
      std::ios_base::openmode which = std::ios_base::in | std::ios_base::out )
    {
      setg(  &myBuffer , &myBuffer , &myBuffer  ) ;
      //slight hack - if mySource has read past the end, won't let us seek
      //and there's no good way to propagate a clear to the original stream
      mySource->clear();
      mySource->seekg(sp);
      std::streampos ret = mySource->tellg();
//      sync();
      return ret;
    };

    /// Returns current source.
    std::istream* GetSource()const
    {
      return mySource;
    };

    ///Changes the source
    void SetSource(std::istream* newsource)
    {
      mySource = newsource;
      setg( &myBuffer , &myBuffer , &myBuffer + 1 ) ;
    }

//    Extractor&   extractor() {return myExtractor;};

  private:
    std::istream*          mySource ;
    Extractor                myExtractor ;
    char                     myBuffer ;
  } ;

//*******************************************************
  template< class Extractor >
  FilteringInputStreambuf< Extractor >::FilteringInputStreambuf(
    std::istream*        source )
    : mySource(source)
  {
    setg( &myBuffer , &myBuffer , &myBuffer ) ;
  }

  //////////////////////////////////////////////////////////
  template< class Extractor >
  int
  FilteringInputStreambuf< Extractor >::underflow()
  {
    int result( EOF ) ;
    if( gptr() < egptr() )
      result = *gptr() ;
    else if ( mySource != nullptr )
    {
      result = myExtractor( *mySource ) ;
      if ( result != EOF )
      {
        if( result < 0 || result > UCHAR_MAX )
          std::cerr << "FilteringInputStreambuf error" << std::endl;
        myBuffer = result ;
        setg( &myBuffer , &myBuffer , &myBuffer + 1 ) ;
      }
    }
    return result ;
  }

  //////////////////////////////////////////////////////
  template< class Extractor >
  int
  FilteringInputStreambuf< Extractor >::sync()
  {
    int result( 0 ) ;
    if ( mySource != nullptr )
    {
      if ( gptr() < egptr() )
      { //have read something but not provided it
        mySource->putback(*gptr() ) ;
        setg( &myBuffer , &myBuffer , &myBuffer ) ;
      }
      result = mySource->sync();
    }
    return result ;
  }

//*********************************************
/// \class LineEndingExtractor lineend.h <openbabel/lineend.h>
  /// \brief Replaces CRLF (DOS) and CR (Mac OS 9) line endings by LF (POSIX)
class OBCONV LineEndingExtractor
{
public:
  int operator()( std::istream& src )
  {
    int ch( src.get() ) ;
    switch (ch)
    {
      case 13: //CR or CRLF
        if(src.peek() == 10)
          src.get(); //CRLF
        //fall through
      case 10: //LF
        return '\n';
        break;
      default:
        return ch;
    }
  }
  void finalize( std::streambuf& ) {}
};

/*! \class FilteringInputStream lineend.h <openbabel/lineend.h>
  \brief A stream interface for FilteringInputStreambuf
  */
template <class Extractor >
class  FilteringInputStream :
    virtual private FilteringInputStreambuf<Extractor>,
    public std::istream
{
public:
    typedef std::istream& istream_reference;
    typedef std::istream istream_type;

    explicit FilteringInputStream(istream_reference istream):
        FilteringInputStreambuf<Extractor>(&istream),std::istream(this) {}
    virtual ~FilteringInputStream() {}

};

} //namespace

#endif //OB_LINEEND_H
//! \file lineend.h
//! \brief Translate line endings automatically (UNIX, Classic Mac, DOS)

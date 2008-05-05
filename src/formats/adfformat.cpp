//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//

// STD
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>


#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include "t41data.h"

using namespace std;
using namespace OpenBabel;


static const double BOHR_TO_ANGSTROM = 0.529177249;

namespace OpenBabel {

class OBT41Format : public OBMoleculeFormat
{
public:
    /// Constructor: register 't41' and "T41" format.
    OBT41Format()
    {
        OBConversion::RegisterFormat( "t41", this );
        OBConversion::RegisterFormat( "T41", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "ADF ASCII t41 format\n"
        "Read only.\n"
        "b no bonds\n"
        "s no multiple bonds\n\n";
    }

    /// Return a specification url, not really a specification since
    /// I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.scm.com/Doc/Doc2006.01/ADF/Analysis/page8.html";
    }

    /// Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

      /// Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return READONEONLY;
    };

    /// Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OBConversion* pConv ) { return 0; }

    /// Read.
    virtual bool ReadMolecule( OBBase* pOb, OBConversion* pConv );

    /// Write: always returns false.
    virtual bool WriteMolecule( OBBase* , OBConversion* )
    {
          return false;
    }

private:
    ///Maps atom name to atomic number.
    int GetAtomicNumber( const string& name ) const
    {
        int iso;
        return etab.GetAtomicNum( name.c_str(), iso );
    }

    ///Utility function that eats all the remaining characters on the current and next line.
    void eol( istream& is ) const { string s; getline( is, s ); getline( is, s ); }

    ///Advance to next tag.
    bool NextTag( istream& is, const std::string& tag ) const
    {
        string buf = "";
        while( is >> buf ) if( buf == tag ) return true;
        return false;
    }

    /// 3x3 Matrix - 3d vector inplace multiply.
    void MatVecMul( const double xColumn[ 3 ],
                    const double yColumn[ 3 ],
                    const double zColumn[ 3 ],
                    double v[ 3 ] )
    {
        double t[ 3 ];
        t[ 0 ] = v[ 0 ];
        t[ 1 ] = v[ 1 ];
        t[ 2 ] = v[ 2 ];
        v[ 0 ] = xColumn[ 0 ] * t[ 0 ] + yColumn[ 0 ] * t[ 1 ] + zColumn[ 0 ] * t[ 2 ];
        v[ 1 ] = xColumn[ 1 ] * t[ 0 ] + yColumn[ 1 ] * t[ 1 ] + zColumn[ 1 ] * t[ 2 ];
        v[ 2 ] = xColumn[ 2 ] * t[ 0 ] + yColumn[ 2 ] * t[ 1 ] + zColumn[ 2 ] * t[ 2 ];
    }


    /// Add grids from SCF
    void AddSCFGrids( istream& is, OBT41Data& t41 ) {}

    ///Inner class used to hold atomic number, coordinate, charge data
    struct AtomData
    {
        int atomicNum;
        double coord[ 3 ];
        double charge;
        AtomData() : atomicNum( 0 ), charge( 0. ) {}
        AtomData( int an ) : atomicNum( an ) {}
        AtomData( int an, const double ac[ 3 ], double c ) : atomicNum( an ), charge( c )
        {
            coord[ 0 ] = ac[ 0 ];
            coord[ 1 ] = ac[ 1 ];
            coord[ 2 ] = ac[ 2 ];
        }
    };

    ///Inner class used to store grid data info before adding it
    ///into an OBMol instance.
    struct T41GridData
    {
        bool valid;
        double startPoint[ 3 ];
        int numPoints[ 3 ];
        double xAxis[ 3 ];
        double yAxis[ 3 ];
        double zAxis[ 3 ];
        int numSymmetries;
        std::vector< std::string > labels;
        bool unrestricted;
        T41GridData() : valid( false ) {}
        operator bool() const { return valid; }
    };

    typedef T41GridData GridData;

    ///Read grid data.
    GridData ReadGridData( istream& is ) const;

    ///Read SCF grids.
    bool ReadSCFGrid( istream& is, OBT41Data& t41Data ) const;
    
    ///Read SCF orbital grids.
    bool ReadSCFOrbitalGrid( istream& is, OBT41Data& t41Data ) const;

    ///Read SumFrag grids.
    bool ReadSumFragGrid( istream& is, OBT41Data& t41Data ) const;

};

//------------------------------------------------------------------------------

namespace
{
    // Global variable used to register Tape41 format.
    OBT41Format t41Format__;
}

//------------------------------------------------------------------------------


//==============================================================================

//------------------------------------------------------------------------------
bool OBT41Format::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    istream& ifs = *pConv->GetInStream();

    GridData gd;
    gd = ReadGridData( ifs );

    OBT41Data* t41Data = 0;
    if( gd )
    {
       // vector< Orbital > orbitals;
       // ReadOrbitals( orbitals, gd.labels[ 0 ] );
       t41Data = new OBT41Data;
       t41Data->SetNumberOfPoints( gd.numPoints[ 0 ], gd.numPoints[ 1 ], gd.numPoints[ 2 ] );
       t41Data->SetAxes( gd.xAxis, gd.yAxis, gd.zAxis );
       t41Data->SetStartPoint( gd.startPoint );
       t41Data->SetUnrestricted( gd.unrestricted );
       t41Data->SetNumSymmetries( gd.numSymmetries );
       streampos current = ifs.tellg();
       while( ReadSCFOrbitalGrid( ifs, *t41Data ) );
       ifs.clear();
       ifs.seekg( current, ios::beg );
       while( ReadSCFGrid( ifs, *t41Data ) );
       ifs.clear();
       ifs.seekg( current, ios::beg );
       while( ReadSumFragGrid( ifs, *t41Data ) );
       pmol->SetData( t41Data );
       ifs.clear();
       ifs.seekg( current, ios::beg );
    }

    string buf;
    // nuuc
    while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
    ifs >> buf; cout << buf << endl;
    if( buf != "nnuc" )
    {
        obErrorLog.ThrowError( __FUNCTION__, "no 'nuuc' after first Geometry tag" );
        return false;
    }
    eol( ifs );
    int numAtoms = -1;
    ifs >> numAtoms; cout << numAtoms << endl;
    buf  = "";

    // labels
    while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
    ifs >> buf; cout << buf << endl;
    if( buf != "labels" )
    {
        obErrorLog.ThrowError( __FUNCTION__, "no 'labels' after second Geometry tag" );
        return false;
    }
    eol( ifs );
    std::vector< AtomData > atoms;
    atoms.reserve( numAtoms );
    for( int i = 0; i != numAtoms; ++i )
    {
        ifs >> buf; cout << buf << endl;
        atoms.push_back( GetAtomicNumber( buf ) );
    }
    if( atoms.size() != numAtoms )
    {
        obErrorLog.ThrowError( __FUNCTION__, "wrong number of atoms" );
        return false;
    }
    //coordinates
    buf = "";
    while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
    ifs >> buf; cout << buf << endl;
    if( buf != "xyznuc" )
    {
        obErrorLog.ThrowError( __FUNCTION__, "no 'xyznuc' after third Geometry tag" );
        return false;
    }
    eol( ifs );
    for( int i = 0; i != numAtoms; ++i )
    {
        ifs >> atoms[ i ].coord[ 0 ] >> atoms[ i ].coord[ 1 ] >> atoms[ i ].coord[ 2 ];
        cout << atoms[ i ].coord[ 0 ] << ' ' << atoms[ i ].coord[ 1 ] << ' ' << atoms[ i ].coord[ 2 ] << endl;
    }
    //charge
    buf = "";
    while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
    ifs >> buf; cout << buf << endl;
    if( buf != "qtch" )
    {
        obErrorLog.ThrowError( __FUNCTION__, "no 'qtch' after fourth Geometry tag" );
        return false;
    }
    eol( ifs );
    for( int i = 0; i != numAtoms; ++i )
    {
        ifs >> atoms[ i ].charge;
    }

    // unit of length
    buf = "";
    while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
    ifs >> buf >> buf >> buf; cout << buf << endl;
    if( buf != "length" )
    {
        obErrorLog.ThrowError( __FUNCTION__, "no 'unit of length' after fifth Geometry tag" );
        return false;
    }
    eol( ifs );
    double scale = 1.0;
    ifs >> scale; 
    /// @todo multply coordinates by axis length;
    for( int i = 0; i != numAtoms; ++i )
    {

        atoms[ i ].coord[ 0 ] *= BOHR_TO_ANGSTROM;
        atoms[ i ].coord[ 1 ] *= BOHR_TO_ANGSTROM;
        atoms[ i ].coord[ 2 ] *= BOHR_TO_ANGSTROM;
//        atoms[ i ].coord[ 0 ] *= scale;
//        atoms[ i ].coord[ 1 ] *= scale;
//        atoms[ i ].coord[ 2 ] *= scale;
//        MatVecMul( gd.xAxis, gd.yAxis, gd.zAxis, atoms[ i ].coord );
//        atoms[ i ].coord[ 0 ] += gd.startPoint[ 0 ];
//        atoms[ i ].coord[ 1 ] += gd.startPoint[ 1 ];
//        atoms[ i ].coord[ 2 ] += gd.startPoint[ 2 ];
    }

    // build OB molecule

    pmol->BeginModify();

    pmol->SetDimension( 3 );

    pmol->ReserveAtoms( numAtoms );

    for( int i = 0; i < numAtoms; ++i )
    {
        OBAtom *atom = pmol->NewAtom();
        atom->SetAtomicNum( atoms[ i ].atomicNum );
        atom->SetVector( atoms[ i ].coord[ 0 ],
                         atoms[ i ].coord[ 1 ],
                         atoms[ i ].coord[ 2 ] );
        atom->SetPartialCharge( atoms[ i ].charge );
    }

    if( !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) pmol->ConnectTheDots();
    if (!pConv->IsOption( "s", OBConversion::INOPTIONS )
        && !pConv->IsOption( "b", OBConversion::INOPTIONS ) )
    {
        pmol->PerceiveBondOrders();
    }
    pmol->EndModify();

    return true;
}

//------------------------------------------------------------------------------
OBT41Format::GridData OBT41Format::ReadGridData( istream& is ) const
{
    GridData gd;
    string buf;
    // Start_point
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "Start_point" ) return gd;
    eol( is );
    is >> gd.startPoint[ 0 ] >> gd.startPoint[ 1 ] >> gd.startPoint[ 2 ];

    gd.startPoint[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.startPoint[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.startPoint[ 2 ] *= BOHR_TO_ANGSTROM;

    // nr of points x
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "x" ) return gd;
    eol( is );
    is >> gd.numPoints[ 0 ];
    // nr of points y
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "y" ) return gd;
    eol( is );
    is >> gd.numPoints[ 1 ];
    // nr of points z
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "z" ) return gd;
    eol( is );
    is >> gd.numPoints[ 2 ];
    // total nr of points
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "points" ) return gd;
    eol( is );
    int n = 0;
    is >> n;
    if( gd.numPoints[ 0 ] * gd.numPoints[ 1 ] * gd.numPoints[ 2 ] != n ) return gd;
    //x-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "x-vector" ) return gd;
    eol( is );
    is >> gd.xAxis[ 0 ] >> gd.xAxis[ 1 ] >> gd.xAxis[ 2 ];

    gd.xAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.xAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.xAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //y-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "y-vector" ) return gd;
    eol( is );
    is >> gd.yAxis[ 0 ] >> gd.yAxis[ 1 ] >> gd.yAxis[ 2 ];

    gd.yAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.yAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.yAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //z-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "z-vector" ) return gd;
    eol( is );
    is >> gd.zAxis[ 0 ] >> gd.zAxis[ 1 ] >> gd.zAxis[ 2 ];

    gd.zAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.zAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.zAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //nr of symmetries
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf;
    if( buf != "symmetries" ) return gd;
    eol( is );
    is >> gd.numSymmetries;
    //labels ///@warning only one label supported
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "labels" ) return gd;
    eol( is );
    is >> buf; gd.labels.push_back( buf );
    //unrestricted
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "unrestricted" ) return gd;
    eol( is );
    char c;
    is >> c;
    gd.unrestricted = ( c == 'T' );

    gd.valid = true;
    return gd;
}

//------------------------------------------------------------------------------
inline bool IsNum( const string& s )
{
    bool isnum = true;
    for( int i = 0; i != s.size(); ++i )
    {
        if( !isdigit( s[ i ] ) )
        {
            isnum = false;
            break;
        }
    }
    return isnum;
}

bool OBT41Format::ReadSCFOrbitalGrid( istream& is, OBT41Data& t41Data ) const
{
    //find next tag starting with 'SCF'
    //if tag starts with SCF_ check next line
    //  advance to next SCF_ tag until tag on next line is a number
    //  read number = orbital id
    //  skip line
    //  read grid data
    //  add values to t41Data: t41Data.SetValues( tag + ' ' + number, double vector );
    if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf.find( "SCF", 0 ) == 0 && buf.size() > 3 ) break;
    if( !is ) return false;
    // SCF_
    const string scf = buf;
    buf = "";
    is >> buf;
    if( !IsNum( buf ) ) // number ?
    {
        while( is >> buf ) // not a number keep on reading
        {
            if( buf == scf ) // SCF_X_X
            {
                is >> buf; // read tag on next line
                if( IsNum( buf ) ) break; // break if (orbital) number
            }
        }
    }
    if( !is ) return false; // eof -> return
    // if we get here it means we are past the orbital number so read
    // read grid values
    const string label = scf + ' ' + buf; cout << label << endl;
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    eol( is );
    if( !is ) return false;
    for( int i = 0; i != numPoints; ++i )
    {
        is >> grid[ i ];
    }
    t41Data.SetValues( label, grid );
    return true;
}

//------------------------------------------------------------------------------
bool OBT41Format::ReadSCFGrid( istream& is, OBT41Data& t41Data ) const
{
	if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf.find( "SCF", 0 ) == 0 && buf.size() == 3 ) break;
    if( !is ) return false;
    // if tag = SCF read next line then skip line and read grid data
    const string scf = buf; // SCF
    is >> buf; // tag on line after SCF
    const string label = scf + ' ' + buf; cout << label << endl;
    eol( is );
    if( !is ) return false;
    // read grid data
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    int i = 0;
    for( ; i != numPoints; ++i ) is >> grid[ i ];
    t41Data.SetValues( label, grid );
    return true;   
}


//------------------------------------------------------------------------------
bool OBT41Format::ReadSumFragGrid( istream& is, OBT41Data& t41Data ) const
{
    if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf == "SumFrag" ) break; // look for SumFrag string
    if( !is ) return false; // not found -> return
    const string sumfrag = buf; // found read next line then skip one line and read data
    is >> buf;
    const string label = sumfrag + ' ' + buf; cout << label << endl;
    eol( is );
    if( !is ) return false;
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    int i = 0;
    for( ; i != numPoints; ++i ) is >> grid[ i ];
    t41Data.SetValues( label, grid );
    return true;
}

} // end namespace OpenBabel

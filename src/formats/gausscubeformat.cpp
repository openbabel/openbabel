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
// reference: http://www.gaussian.com/g_ur/u_cubegen.htm

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/griddata.h>

using namespace std;
using namespace OpenBabel;


/// @warning seems that every position in cube files has always to be
/// multiplied by this constant even if according to:
/// reference: http://www.gaussian.com/g_ur/u_cubegen.htm
/// <Num points along first axis> > 0 means that the values
/// are already in Angstroms.
static const double BOHR_TO_ANGSTROM = 0.529177249;

/// @todo remove AtomPostion and GaussianCube structs and implement everything with OBGridData/OBMol without
/// using any other data type.
/// @todo add error checking code.
namespace
{
    struct AtomPosition
    {
        int atomicNumber;
        double value; /// @todo store this value into OBMol somehow, consider using per-atom
        double position[ 3 ];
    };

    struct GaussianCube
    {
        typedef enum { BOHR, ANGSTROM } Unit;
        std::string firstLine;
        std::string secondLine;
        int numberOfAtoms;
        double origin[ 3 ];
        int numPoints[ 3 ];
        double xAxis[ 3 ];
        double yAxis[ 3 ];
        double zAxis[ 3 ];
        std::vector< AtomPosition > atomPositions;
        std::vector< double > values;
        Unit unit;
    };

    bool ReadGaussianCube( GaussianCube& gc, istream& in )
    {
        try
        {
            static const int MAX_LINE_SIZE = 512;
            char buffer[ MAX_LINE_SIZE ];
            // read title
            in.getline( buffer, MAX_LINE_SIZE );
            gc.firstLine = string( buffer );
            in.getline( buffer, MAX_LINE_SIZE );
            gc.secondLine = string( buffer );
            // number of atoms and center
            in >> gc.numberOfAtoms >> gc.origin[ 0 ] >> gc.origin[ 1 ] >> gc.origin[ 2 ];

            gc.origin[ 0 ] *= BOHR_TO_ANGSTROM;
            gc.origin[ 1 ] *= BOHR_TO_ANGSTROM;
            gc.origin[ 2 ] *= BOHR_TO_ANGSTROM;

            bool negativeNumAtoms = false;
            if( gc.numberOfAtoms < 0 )
            {
                 gc.numberOfAtoms = -gc.numberOfAtoms;
                 negativeNumAtoms = true; // should mean that there is some information
                                          // between atom positions and grid values.
            }
            // point number and axis
            in >> gc.numPoints[ 0 ] >> gc.xAxis[ 0 ] >> gc.xAxis[ 1 ] >> gc.xAxis[ 2 ];
            in >> gc.numPoints[ 1 ] >> gc.yAxis[ 0 ] >> gc.yAxis[ 1 ] >> gc.yAxis[ 2 ];
            in >> gc.numPoints[ 2 ] >> gc.zAxis[ 0 ] >> gc.zAxis[ 1 ] >> gc.zAxis[ 2 ];

            gc.xAxis[ 0 ] *= BOHR_TO_ANGSTROM; gc.xAxis[ 1 ] *= BOHR_TO_ANGSTROM; gc.xAxis[ 2 ] *= BOHR_TO_ANGSTROM;
            gc.yAxis[ 0 ] *= BOHR_TO_ANGSTROM; gc.yAxis[ 1 ] *= BOHR_TO_ANGSTROM; gc.yAxis[ 2 ] *= BOHR_TO_ANGSTROM;
            gc.zAxis[ 0 ] *= BOHR_TO_ANGSTROM; gc.zAxis[ 1 ] *= BOHR_TO_ANGSTROM; gc.zAxis[ 2 ] *= BOHR_TO_ANGSTROM;

            // set unit
            gc.unit = gc.numPoints[ 0 ] < 0 ? GaussianCube::BOHR : GaussianCube::ANGSTROM;
            gc.numPoints[ 0 ] = abs( gc.numPoints[ 0 ] );
            gc.atomPositions.reserve( gc.numberOfAtoms  );
            // read atoms
            for( int i = 0; i < gc.numberOfAtoms; ++i )
            {
                AtomPosition ap;
                in >> ap.atomicNumber >> ap.value >> ap.position[ 0 ]
                   >> ap.position[ 1 ] >> ap.position[ 2 ];

                ap.position[ 0 ] *= BOHR_TO_ANGSTROM;
                ap.position[ 1 ] *= BOHR_TO_ANGSTROM;
                ap.position[ 2 ] *= BOHR_TO_ANGSTROM;
                gc.atomPositions.push_back( ap );
            }

            // if number of atoms is negative means that there is some data between
            // atom positions and grid data
            if( negativeNumAtoms )
            {
                int n = 0;
                in >> n;
                int dummy;
                for( int j = 0; j < n; ++j ) in >> dummy;
            }

            // read values
            gc.values.reserve( gc.numPoints[ 0 ] *
                               gc.numPoints[ 1 ] *
                               gc.numPoints[ 2 ] );
            while( in )
            {
                double v;
                in >> v;
                gc.values.push_back( v );
            }
        }
        catch( const exception& )
        {
            return false;
        }

        return true;
    }
}

//==============================================================================
/// Class to read Gaussian cube files.
/// Atoms are stored into an instance of OBMol and grid data into an instance
/// of OBGridData, a sublass of OBGenericData.
/// To retrieve the grid data from an OBMol instance:
/// @code
/// OBMol* mol;
/// ...
/// if( mol->HasData( "GridData" ) )
/// {
///    const OBGridData* gd = dynamic_cast< const OBGridData* >( mol->GetData( "GridData" ) );
///    assert( gd );
///    ...
/// }
/// @endcode
/// @note Instead of "GridData" OBGenericDataType::CustomData0 can be used instead.
///
/// Sample code to print the data contained into an OBGridData instance.
/// @code
/// void PrintGaussianCubeData( const OBGenericData* d )
///{
///	const OBGridData* gd = dynamic_cast< const OBGridData* >( d );
///	if( !gd )
///	{
///		cerr << "WRONG TYPE" << endl;
///		return;
///	}
///
///
///	cout << "UNIT:            ";
///	if( gd->GetUnit() == OBGridData::BOHR ) cout << "BOHR";
///	else cout << "ANGSTROM";
///	cout << '\n';
///
///	cout << "ORIGIN:          " << gd->GetOrigin()[ 0 ]   << ' '
///								<< gd->GetOrigin()[ 1 ]   << ' '
///								<< gd->GetOrigin()[ 2 ]   << '\n';
///	int nx, ny, nz;
///	gd->GetNumberOfPoints( nx, ny, nz );
///	cout << "NUM POINTS:      "	<< nx << ' '
///								<< ny << ' '
///								<< nz << '\n';
///	double xAxis[ 3 ], yAxis[ 3 ], zAxis[ 3 ];
///	gd->GetAxes( xAxis, yAxis, zAxis );
///	cout << "X AXIS:          " << xAxis[ 0 ] << ' '
///								<< xAxis[ 1 ] << ' '
///								<< xAxis[ 2 ] << '\n';
///	cout << "Y AXIS:          " << yAxis[ 0 ] << ' '
///								<< yAxis[ 1 ] << ' '
///								<< yAxis[ 2 ] << '\n';
///	cout << "Z AXIS:          " << zAxis[ 0 ] << ' '
///								<< zAxis[ 1 ] << ' '
///								<< zAxis[ 2 ] << '\n';
///	copy( gd->GetValues().begin(), gd->GetValues().end(), ostream_iterator< double >( cout, "\n" ) );
///
///}
///
/// @endcode
/// @todo add support for optional premultiplication of atom positions by BOHR_TO_ANGSTROM constant.
class OBGaussianCubeFormat : public OpenBabel::OBMoleculeFormat
{
public:
    /// Constructor: register 'cube' and "CUBE" format.
    OBGaussianCubeFormat()
    {
        OpenBabel::OBConversion::RegisterFormat( "cube", this );
        OpenBabel::OBConversion::RegisterFormat( "CUBE", this );
        OpenBabel::OBConversion::RegisterFormat( "cub", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "Gaussian cube format\n"
        "Read only.\n"
        "b no bonds\n"
        "s no multiple bonds\n\n";
    }

    /// Return a specification url, not really a specification since
    /// I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.gaussian.com/g_ur/u_cubegen.htm";
    }

    /// Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

      /// Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return READONEONLY;
    };

    /// Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    /// Read.
    virtual bool ReadMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );

    /// Write: always returns false.
    virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* )
    {
        return false;
    }
};

//------------------------------------------------------------------------------

namespace
{
    // Global variable used to register Gaussian cube format.
    OBGaussianCubeFormat gaussianCubeFormat__;
}

//------------------------------------------------------------------------------


//==============================================================================

//------------------------------------------------------------------------------
bool OBGaussianCubeFormat::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    istream& ifs = *pConv->GetInStream();

    pmol->BeginModify();

    /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

    GaussianCube gc;
    if( !ReadGaussianCube( gc, ifs ) )
    {
        obErrorLog.ThrowError( __FUNCTION__, "Problems reading a Gaussian cube file.", obWarning);
        return false;
    }

    pmol->SetDimension( 3 );
    pmol->ReserveAtoms( gc.numberOfAtoms );
    pmol->SetTitle( gc.firstLine );
    for( int i = 0; i < gc.numberOfAtoms; ++i )
    {
        OBAtom *atom = pmol->NewAtom();
        atom->SetAtomicNum( gc.atomPositions[ i ].atomicNumber );
        atom->SetVector( gc.atomPositions[ i ].position[ 0 ],
                         gc.atomPositions[ i ].position[ 1 ],
                         gc.atomPositions[ i ].position[ 2 ] );
        // use OBPairData to store value read from file before atom position.
        OBPairData* v = new OBPairData;
        v->SetAttribute( "GaussianCubeValue" );
        ostringstream ss;
        ss << gc.atomPositions[ i ].value;
        const string val = ss.str();
        char* strBuffer = new char[ val.size() + 1 ];
        strcpy( strBuffer, val.c_str() );
        strBuffer[ val.size() ] = 0;
        v->SetValue( ss.str().c_str() );
        atom->SetData( v );
    }

    OBGridData* gd = new OBGridData;

    gd->SetNumberOfPoints( gc.numPoints[ 0 ],
                           gc.numPoints[ 1 ],
                           gc.numPoints[ 2 ] );
    gd->SetLimits( gc.origin, gc.xAxis, gc.yAxis, gc.zAxis );
    gd->SetUnit( gc.unit == GaussianCube::BOHR ? OBGridData::BOHR : OBGridData::ANGSTROM );
    gd->SetValues( gc.values );

    cerr << "min " <<  *std::min_element( gc.values.begin(), gc.values.end() ) << endl;
    cerr << "max " <<  *std::max_element( gc.values.begin(), gc.values.end() ) << endl;

    gd->SetOrigin(fileformatInput); // i.e., is this data from a file or determined by Open Babel

    if( !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) pmol->ConnectTheDots();
    if (!pConv->IsOption( "s", OBConversion::INOPTIONS )
        && !pConv->IsOption( "b", OBConversion::INOPTIONS ) )
    {
        pmol->PerceiveBondOrders();
    }
    pmol->EndModify();

    pmol->SetData( gd );

    return true;
}

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
// reference: http://www.cmbi.ru.nl/molden/molden_format.html

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

using namespace std;

namespace OpenBabel
{

/// Molden input reader: reads atoms from [Atoms] section of Molden input file.
class OBMoldenFormat : public OpenBabel::OBMoleculeFormat
{
public:
    /// Constructor: register 'molden' format.
    OBMoldenFormat()
    {
        OBConversion::RegisterFormat( "molden", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "Molden input format\n"
        "ReadOnly.\n"
        "Read Options e.g. -as\n"
        "  b no bonds\n"
        "  s no multiple bonds\n\n";
    }

    /// Return a specification url, not really a specification since
    /// I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.cmbi.ru.nl/molden/molden_format.html";
    }

    /// Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

      /// Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return READONEONLY | NOTWRITABLE  ;
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
    // Global variable used to register Molden format.
    OBMoldenFormat moldenFormat__;
}

//------------------------------------------------------------------------------


//==============================================================================

//------------------------------------------------------------------------------
bool OBMoldenFormat::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    istream& ifs = *pConv->GetInStream();

    pmol->BeginModify();
    pmol->SetDimension( 3 );
    string lineBuffer;
    getline( ifs, lineBuffer );

    while( ifs && lineBuffer.find( "[Atoms]" ) == string::npos )
    {
        getline( ifs, lineBuffer );
    }

    if( !ifs ) return false;

    //[Atoms] AU OR Angs
    double factor = 1.; // Angstrom
    if( lineBuffer.find( "AU" ) != string::npos ) factor = 0.529177249; // Bohr

    while( ifs )
    {
        getline( ifs, lineBuffer );
        if( lineBuffer.size() == 0 ) continue;
        if( lineBuffer.find( '[' ) != string::npos ) break;
        istringstream is( lineBuffer );
        string atomName;
        int atomId;
        int atomicNumber;
        double x, y, z;
        is >> atomName >> atomId >> atomicNumber >> x >> y >> z;
        OBAtom* atom = pmol->NewAtom();
        if( !atom ) break;
        atom->SetAtomicNum( atomicNumber );
        atom->SetVector( x * factor, y * factor, z * factor );
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

}

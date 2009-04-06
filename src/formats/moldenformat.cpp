//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)
// Some portions Copyright (C) 2009 Michael Banck
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

    //Vibrational data
    std::vector< std::vector< vector3 > > Lx;
    std::vector<double> Frequencies, Intensities;

    pmol->BeginModify();
    pmol->SetDimension( 3 );
    string lineBuffer;

    while( getline( ifs, lineBuffer ) )
      {
        if( lineBuffer.find( "[Atoms]" ) != string::npos || 
            lineBuffer.find( "[ATOMS]" ) != string::npos ) {
          double factor = 1.; // Angstrom
          if( lineBuffer.find( "AU" ) != string::npos ) factor = 0.529177249; // Bohr
          getline( ifs, lineBuffer );
          while( lineBuffer.find( "[") == string::npos )
            {
              if( lineBuffer == "" ) continue;
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
              getline( ifs, lineBuffer );
            }
	  continue;
        } // "[Atoms]" || "[ATOMS]"
        if( lineBuffer.find( "[FREQ]" ) != string::npos ) {        
          while( getline( ifs, lineBuffer ) )
            {
              if( lineBuffer == "" ) continue;
              if( lineBuffer.find( "[" ) != string::npos ) break;
              istringstream is( lineBuffer );
              double freq;
              is >> freq;
              Frequencies.push_back( freq );
            }
          continue;
        } // "[FREQ]"
        if( lineBuffer.find( "[INT]" ) != string::npos ) {        
          while( getline( ifs, lineBuffer ) )
            {
              if( lineBuffer == "" ) continue;
              if( lineBuffer.find( "[" ) != string::npos ) break;
              istringstream is( lineBuffer );
              double intens;
              is >> intens;
              Intensities.push_back( intens );
            }
          continue;
        } // "[INT]"
        if( lineBuffer.find( "[FR-NORM-COORD]" ) != string::npos ) {
          getline( ifs, lineBuffer );
          while( ifs && lineBuffer.find( "Vibration") != string::npos ) 
            {
              vector<vector3> vib;
              getline( ifs, lineBuffer );
              while( ifs && lineBuffer.find( "Vibration") == string::npos ) 
                {
                  istringstream is( lineBuffer );
                  double x, y, z;
                  is >> x >> y >> z;
                  vib.push_back( vector3( x, y, z ) );
                  getline( ifs, lineBuffer );
                }
              Lx.push_back( vib );  
           } // while
        } // "[FR-NORM-COORD]"
      } // while

    if ( pmol->NumAtoms() == 0 ) {
      pmol->EndModify();
      return false;
    }

    // Attach vibrational data, if there is any, to molecule
    if(Frequencies.size()>0)
    {
      for (int i=0; i < Frequencies.size(); i++) {
        // Set intensities to zero if none were specified
        if (Intensities.size() < Frequencies.size() ) Intensities.push_back( 0.0 );
        if (fabs(Frequencies[i]) < 10.) {
          // skip translational and rotational modes
          Frequencies.erase( Frequencies.begin() + i );
          Intensities.erase( Intensities.begin() + i );
          Lx.erase( Lx.begin() + i );  
          i--;  // compensate for the vibration which just got cut out
        }
      }
      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(Lx, Frequencies, Intensities);
      pmol->SetData(vd);
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

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

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

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
        OBConversion::RegisterFormat( "mold", this );
        OBConversion::RegisterFormat( "molf", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "Molden format\n"
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

      /// Return read/write flag.
    virtual unsigned int Flags()
    {
        return READONEONLY | WRITEONEONLY ;
    };

    /// Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    /// Read.
    virtual bool ReadMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );

    /// Write.
    virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
};

//------------------------------------------------------------------------------

// Global variable used to register Molden format.
OBMoldenFormat moldenFormat__;

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
          if( lineBuffer.find( "AU" ) != string::npos ) factor = BOHR_TO_ANGSTROM; // Bohr
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
        } // "[INT]"
        if( lineBuffer.find( "[FR-COORD]" ) != string::npos ) {
          if (pmol->NumAtoms() == 0) {
            // No atoms yet, probably there is no [ATOMS] section
            // in the file.
            while ( getline( ifs, lineBuffer ) )
              {
                if( lineBuffer == "" ) continue;
                if( lineBuffer.find( "[" ) != string::npos ) break;
                string atomName;
                double x, y, z;
                istringstream is( lineBuffer );
                is >> atomName >> x >> y >> z;
                OBAtom* atom = pmol->NewAtom();
                if( !atom ) break;
                atom->SetAtomicNum( etab.GetAtomicNum(atomName.c_str()));
                // Vibrational equilibrium geometry is mandated to be
                // in Bohr.
                atom->SetVector( x * BOHR_TO_ANGSTROM,
                                 y * BOHR_TO_ANGSTROM,
                                 z * BOHR_TO_ANGSTROM);
             }
           }
         } // "[FR-COORD]"
        if( lineBuffer.find( "[FR-NORM-COORD]" ) != string::npos ) {
          getline( ifs, lineBuffer );
          while( ifs && lineBuffer.find( "ibration") != string::npos )
            {
              vector<vector3> vib;
              getline( ifs, lineBuffer );
              while( ifs && lineBuffer.find( "ibration") == string::npos )
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
        if (fabs(Frequencies[i]) < 10.) {
          // skip translational and rotational modes
          Frequencies.erase( Frequencies.begin() + i );
          if (Intensities.size() > i) Intensities.erase( Intensities.begin() + i );
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

bool OBMoldenFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    int i = 1;

    ofs << "[Molden Format]" << endl;
    ofs << "[Atoms] Angs" << endl;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%2s%6d%3d%13.6f%13.6f%13.6f\n",
                etab.GetSymbol(atom->GetAtomicNum()),
		i++,
                atom->GetAtomicNum(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    OBVibrationData *vib = (OBVibrationData *) mol.GetData(OBGenericDataType::VibrationData);
    if (vib && vib->GetNumberOfFrequencies() > 0) {
      ofs << "[FREQ]" << endl;
      vector<double> frequencies = vib->GetFrequencies();
      vector<double> intensities = vib->GetIntensities();
      for (int i=0; i < vib->GetNumberOfFrequencies(); i++) {
	snprintf(buffer, BUFF_SIZE, "%10.4f\n", frequencies[i]);
        ofs << buffer;
      }
      if (intensities.size() > 0) {
        ofs << "[INT]" << endl;
	for (int i=0; i < vib->GetNumberOfFrequencies(); i++) {
	  snprintf(buffer, BUFF_SIZE, "%10.4f\n", intensities[i]);
	  ofs << buffer;
        }
      }
      ofs << "[FR-COORD]" << endl;
      FOR_ATOMS_OF_MOL(atom, mol)
        {
          snprintf(buffer, BUFF_SIZE, "%2s%13.6f%13.6f%13.6f\n",
                  etab.GetSymbol(atom->GetAtomicNum()),
                  atom->GetX()*ANGSTROM_TO_BOHR,
                  atom->GetY()*ANGSTROM_TO_BOHR,
                  atom->GetZ()*ANGSTROM_TO_BOHR);
          ofs << buffer;
        }
      ofs << "[FR-NORM-COORD]" << endl;
      for (int mode=0; mode < vib->GetNumberOfFrequencies(); mode++) {
	snprintf(buffer, BUFF_SIZE, "vibration%6d\n", mode+1);
	ofs << buffer;
        vector<vector3> lx = vib->GetLx()[mode];
	for (int i=0; i < mol.NumAtoms(); i++) {
	  vector3 disp = lx[i];
	  snprintf(buffer, BUFF_SIZE, "%12.6f%13.6f%13.6f\n",
		  disp[0], disp[1], disp[2]);
	  ofs << buffer;
	}
      }
    } // vib
    return(true);
}

}

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


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>


#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>

using namespace std;

namespace OpenBabel
{


//==============================================================================
/// Class to output a molecule in XYZR MSMS input format for further computation
/// of Connolly surface.
/// Michel Sanner page with info on MSMS:
/// http://www.scripps.edu/~sanner/
class OBMSMSFormat : public OpenBabel::OBMoleculeFormat
{
public:
    /// Constructor: register 'msms' and "MSMS" format.
    OBMSMSFormat()
    {
        OpenBabel::OBConversion::RegisterFormat( "msms", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "M.F. Sanner's MSMS input format\n"
        "Generates input to the MSMS (Michael Sanner Molecular Surface) program to compute solvent surfaces.\n\n"
        "Write Options, e.g. -xa\n"
        "  a output atom names\n";
    }

    /// Return a specification url, not really a specification since
    /// I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.scripps.edu/~sanner";
    }

    /// Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

      /// Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return WRITEONEONLY | NOTREADABLE;
    };

    /// Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    /// Read: always return false.
    virtual bool ReadMolecule( OpenBabel::OBBase*, OpenBabel::OBConversion* )
    {
        return false;
    }

    /// Write.
    virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
};

//------------------------------------------------------------------------------

// Global variable used to register MSMS format.
OBMSMSFormat msmsFormat__;

//------------------------------------------------------------------------------


//==============================================================================

//------------------------------------------------------------------------------
bool OBMSMSFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if( pmol == 0 ) return false;

    ostream& os = *pConv->GetOutStream();

    const bool atomNames = pConv->IsOption( "a", OBConversion::OUTOPTIONS )!=NULL;

    // write header ?

    // iterate through atoms and write <atom x> <atom y> <atom z> <atom radius>
    // and optionally <atomic number> in case atomNames == true

    FOR_ATOMS_OF_MOL( a, *pmol )
    {
        const double* c = a->GetCoordinate();
        os << c[ 0 ] << '\t' << c[ 1 ] << '\t' << c[ 2 ] << '\t' <<
        OBElements::GetVdwRad( a->GetAtomicNum() );
        if( atomNames ) os << '\t' << a->GetAtomicNum();
        os << '\n';
    }
    os.flush();
    return true;
}

}

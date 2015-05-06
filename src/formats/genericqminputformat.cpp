//
// Generic QM Input file writer for Open Babel
// Copyright (C) 2015 M J Harvey, 
// Acellera Ltd 
// m.j.harvey ( at ) acellera.com
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

#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <math.h>

#include "genericqminputformat.h"

using namespace std;


namespace OpenBabel {

	GenericQMInputFormat::GenericQMInputFormat( const char *format_name, const char *type )
	{
		format_type = strdup( format_name );
		OpenBabel::OBConversion::RegisterFormat( type, this );
		OBConversion::RegisterOptionParam("B"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("T"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("E"      , this, 0, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("M"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("O"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("F"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("Q"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("m"      , this, 1, OBConversion::OUTOPTIONS);
		OBConversion::RegisterOptionParam("p"      , this, 1, OBConversion::OUTOPTIONS);

		basis_default = string("6-31g*");
		theory_default= string("rhf");

		basis_sets.insert( std::pair<string,string>(string(basis_default) , string(basis_default )) );
		basis_sets.insert( std::pair<string,string>(string("cc-pvdz"), string ("cc-pvdz" )) );

		theory_sets.insert( std::pair<string,string>(string(theory_default), string (theory_default) ) );
		theory_sets.insert( std::pair<string,string>(string("rohf"), string ("rohf" ) ) );
		theory_sets.insert( std::pair<string,string>(string("uhf"), string ("uhf" ) ) );
		theory_sets.insert( std::pair<string,string>(string("b3lyp"), string ("b3lyp" ) ) );

		mult  = 1;
		espgrid = false;
		mem   = 2;
		ncpus = 1;
		basis = basis_default;
		theory= theory_default;
		opt   = OPT_NONE;
		charge= 0;
		charge_set = false;
		mult_set   = false;
	}

	/// Return description.
	const char* GenericQMInputFormat::Description() 
	{
		stringstream o; 
		o << format_type << " input format \n"
			<< "Write Options:\n"
			<<	 "  T<theory>            Theory level for calculation. Options:   (" << theory_default << ")\n";

		for( std::map<string,string>::iterator i = theory_sets.begin(); i != theory_sets.end(); i++ ) {
			o << "                         " << i->first << "\n";
		}

		o	<<	 "  B<basis-set>         Basis set for calculation. Options:   (" << basis_default << ")\n";

		for( std::map<string,string>::iterator i = basis_sets.begin(); i != basis_sets.end(); i++ ) {
			o << "                         " << i->first << "\n";
		}

		o <<  "  M#                   Multiplicity of the molecule (default to Babel's guess)\n"
			<<  "  Q#                   Charge of the the molecule  (default sum of partial charges)\n"
			<<  "  E                    Calculate field and potential on points from grid.dat\n"
			<<  "  O<type>              Optimise geometry. Options:           (default none)\n"
			<<  "                         none  loose  normal  tight\n"
			<<  "  F<...>               List of bonds, angles and dihedrals to freeze\n"
			<<  "                       Atoms indexed from 1, comma-separated.\n"
			<<  "                       Each term colon-separated.\n"
			<<  "                       Eg freezing one angle, one dihedral: \n"
			<<  "                          -xF 1,2,3:3,5,6,7\n"
			<<  "  m#                   Amount of memory to allocate (GB)     (default 2.0)\n"
			<<  "  p#                   Number of CPUs/threads to allocate    (default 1)\n\n";
		return strdup( o.str().c_str() );
	}

	/// Return read/write flag: read only.
	unsigned int GenericQMInputFormat::Flags()
	{
		return WRITEONEONLY | NOTREADABLE;
	};

	/// Skip to object: used for multi-object file formats.
	int GenericQMInputFormat::SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

	/// Read: always return false.
	bool GenericQMInputFormat::ReadMolecule( OpenBabel::OBBase*, OpenBabel::OBConversion* )
	{
		return false;
	}

	//==============================================================================

	static bool parse_freeze( const char *args_x, int natoms, std::vector< std::vector<int> > &frozen ) {
		char *args = strdup( args_x );
		char *tmp1;
		char *tmp2;
		char *entry;
		bool retval = false;

		entry = strtok_r( args, ":", &tmp1 );
		while( entry != NULL ) {
			char *field = strtok_r( entry, ",", &tmp2 );
			std::vector<int> idx;
			int cnt=0;
			while( cnt<4 && ( field!=NULL ) ) {
				idx.push_back( atoi( field ) );
				if( idx[cnt]<1 || idx[cnt]>natoms ) { goto err; } // invalid index
				cnt++;
				field = strtok_r( NULL, ",", &tmp2 );
			} 
			frozen.push_back( idx );
			if( strtok_r( NULL, ",", &tmp2 ) ) { 
				// one of the entries was too long
				goto err;
			}
			if( idx[0] < 1 ) {
				// one of the entries was too short or invalid index (<1)
				goto err;
			}
			entry = strtok_r( NULL, ":", &tmp1 );
		}

		retval = true;
err:
		free( args );
		return retval;
	}

	//------------------------------------------------------------------------------

	bool GenericQMInputFormat::ParseOptions( OBBase *pOb, OBConversion *pConv ) {
		OBMol* mol = dynamic_cast< OBMol* >(pOb);
		if( mol == 0 ) return false;

		if( pConv->IsOption( "M" , OBConversion::OUTOPTIONS ) ) {
			mult = atoi( pConv->IsOption( "M", OBConversion::OUTOPTIONS) );
			mult_set = true;
			if( mult < 1 ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'M'", obError );
				return false;
			}
		}
		if( pConv->IsOption( "E" , OBConversion::OUTOPTIONS ) ) {
			espgrid = true;
		}
		if( pConv->IsOption( "O" , OBConversion::OUTOPTIONS ) ) {
			const char *o = pConv->IsOption("O", OBConversion::OUTOPTIONS);
			if     ( !strcmp( o, "none"  ) ) { opt = OPT_NONE  ; }
			else if( !strcmp( o, "loose" ) ) { opt = OPT_LOOSE ; }
			else if( !strcmp( o, "normal") ) { opt = OPT_NORMAL; }
			else if( !strcmp( o, "tight" ) ) { opt = OPT_TIGHT ; }
			else if( !strlen( o ) ) { opt = OPT_NORMAL; }
			else {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'O'", obError );
				return false;
			}
		}

		if( pConv->IsOption( "Q" , OBConversion::OUTOPTIONS ) ) {
			charge = atoi( pConv->IsOption( "Q", OBConversion::OUTOPTIONS) );
			charge_set = true;
		} 


		if( pConv->IsOption( "m" , OBConversion::OUTOPTIONS ) ) {
			mem = atoi( pConv->IsOption( "m", OBConversion::OUTOPTIONS) );
			if( mem < 1 ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'm'", obError );
				return false;
			}
		} 
		if( pConv->IsOption( "p" , OBConversion::OUTOPTIONS ) ) {
			ncpus = atoi( pConv->IsOption( "p", OBConversion::OUTOPTIONS) );
			if( ncpus < 1 ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'p'", obError );
				return false;
			}
		} 

		if( pConv->IsOption( "T" , OBConversion::OUTOPTIONS ) ) {
			theory = string( pConv->IsOption("T" , OBConversion::OUTOPTIONS ) );
			if( ! theory_sets.count( theory ) ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'T'", obError );
				return false;
			}
		}


		if( pConv->IsOption( "B" , OBConversion::OUTOPTIONS ) ) {
			basis = string( pConv->IsOption("B" , OBConversion::OUTOPTIONS ) );
			if( ! basis_sets.count( basis ) ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'B'", obError );
				return false;
			}
		}

		if( pConv->IsOption( "F" , OBConversion::OUTOPTIONS ) ) {
			if( ! parse_freeze( pConv->IsOption( "F" , OBConversion::OUTOPTIONS ), mol->NumAtoms(), frozen ) ) {
				obErrorLog.ThrowError(__FUNCTION__, "Invalid value for argument 'F'", obError );
				return false;
			}
		}

		return true;
	}

}

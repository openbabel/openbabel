//
// PSI4GQMInput Input file writer for Open Babel
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

#include "genericqminputformat.h"

#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <math.h>

using namespace std;

namespace OpenBabel
{


	class PSI4GQMInputFormat : public GenericQMInputFormat
	{
		public:
			PSI4GQMInputFormat() : GenericQMInputFormat( "PSI4", "psi4in" ) {
			}

			virtual const char* SpecificationURL() {

				return "http://sirius.chem.vt.edu/psi4manual/latest/index.html";
			}

			virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
	};


	// Global variable used to register PSI4GQMInput format.
	PSI4GQMInputFormat thePSI4GenericQMInputFormat;



	//------------------------------------------------------------------------------
	bool PSI4GQMInputFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
	{
		char line[100];
		if( ! ParseOptions( pOb, pConv ) ) {
			return false;
		}

		OBMol* mol = dynamic_cast< OBMol* >(pOb);
		if( mol == 0 ) return false;

		ostream& os = *pConv->GetOutStream();

		bool is_hf = ( theory.find("hf") != string::npos );

		double q = 0.;

		FOR_ATOMS_OF_MOL(atom, mol) {
			q += atom->GetPartialCharge();
		}

		int qi = (int) round(q);

		if( charge_set ) { qi = charge; }
    if( !mult_set )  { mult = mol->GetTotalSpinMultiplicity(); }


		os << std::setprecision(10);

		if( is_hf ) {
			os << "set reference "<<theory << endl;
		}
		os << "set basis " << basis << endl << endl;



		os << "try:" <<endl;
		os << "\tmem = int(os.getenv(\"MEM\") ) * 1.e9" << endl;
		os << "except:" << endl;
		os << "\tmemory "<< mem << " gb" << endl << endl;

		os << "try:" <<endl;
		os << "\tset_num_threads( int(os.getenv(\"NCPUS\") ) )" << endl ;
		os << "except:" << endl;
		os << "\tset_num_threads(" << ncpus << ")" << endl << endl;

		os << "molecule MOL {" << endl; 
		os << "\t" << qi <<" " << mult << endl;


		FOR_ATOMS_OF_MOL(a, mol) {
			double *c = a->GetCoordinate();
			sprintf( line, "\t%s\t%f\t%f\t%f", etab.GetSymbol(a->GetAtomicNum()) , c[0], c[1], c[2] );
			os << line << endl;
		}
		os << "\tsymmetry c1" << endl;
		os << "}" << endl;

		switch( opt ) {
			case OPT_LOOSE:
				os << "set optking { g_convergence = \"GAU_LOOSE\" }" << endl;
				break;
			case OPT_TIGHT:
				os << "set optking { g_convergence = \"GAU_TIGHT\" }" << endl;
				break;
		}

		if( frozen.size() ) {
			bool first= true;
			os << "set optking { frozen_distance = \" ";
			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {
				if( i->size() == 2 ) {
					if( !first ) { os << " , "; }
					os << (*i)[0] << " " << (*i)[1] ;
					first = false;
				}
			}
			os << " \" }" << endl;

			first=true;
			os << "set optking { frozen_bend = \" ";
			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {
				if( i->size() == 3 ) {
					if( !first ) { os << " , "; }
					os << (*i)[0] << " " << (*i)[1]  << " " << (*i)[2];
					first = false;
				}
			}
			os << " \" }" << endl;

			first=true;
			os << "set optking { frozen_dihedral = \" ";
			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {
				if( i->size() == 4 ) {
					if( !first ) { os << " , "; }
					os << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << " " << (*i)[3];
					first = false;
				}
			}
			os << " \" }" << endl;
		}


		string type = "scf";
		if( !is_hf ) {
			type = theory;
		}
		if( opt == OPT_NONE ) {
			os << "ee = energy('" << type << "')" << endl;
		}
		else {
			os << "ee = optimize('" << type <<"')" << endl;
		}

		// Output coords in XYZ
		os <<  endl << "f = open( 'psi4.xyz', 'w' ) " << endl;
		os << "f.write( \"" << mol->NumAtoms() << " \" )" << endl;
		// nonstandard -- smuggle the QM energy on the first line
		os << "f.write( str(ee) )" << endl;
		os << "f.write( \"\\n\" )" << endl;
		os << "f.write( MOL.save_string_xyz() )" << endl << endl;
		os << "f.close()" << endl << endl;

		if( espgrid ) {
			os << "oeprop( \"GRID_ESP\", \"GRID_FIELD\" )" << endl;
		}
		return true;
	}
}

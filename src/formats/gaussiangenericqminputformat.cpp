//
// GaussianGQMInput Input file writer for Open Babel
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


	class GaussianGQMInputFormat : public GenericQMInputFormat
	{
		public:
			GaussianGQMInputFormat() : GenericQMInputFormat( "Gaussian", "gaussianin" ) {
			}

			virtual const char* SpecificationURL() {

				return "http://www.gaussian.com/g_tech/g09ur.htm";
			}

			virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
	};


	// Global variable used to register GaussianGQMInput format.
	GaussianGQMInputFormat theGaussianGenericQMInputFormat;



	//------------------------------------------------------------------------------
	bool GaussianGQMInputFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
	{
		char line[80];
		if( ! ParseOptions( pOb, pConv ) ) {
			return false;
		}

		OBMol* mol = dynamic_cast< OBMol* >(pOb);
		if( mol == 0 ) return false;

		ostream& os = *pConv->GetOutStream();

		double q = 0.;

		FOR_ATOMS_OF_MOL(atom, mol) {
			q += atom->GetPartialCharge();
		}

		int qi = (int) round(q);

		os << std::setprecision(10);
		os << "%nprocshared="  << ncpus << endl;
		os << "%mem="  << mem << "GB" << endl;
		os << "%chk = chk" << endl;
		os << "%nosave" << endl;
		os << endl;

		// the route

		os << "# " << theory << "/" << basis  << endl; 
		switch( opt ) {
			case OPT_NONE:
				if( frozen.size() ) {
					os << "# opt(modredundant)" << endl;
				}
				break;
			case OPT_NORMAL: 
				os << "# opt" ; 
				if( frozen.size() ) {
					os << "(modredundant)";
				}
				os << endl;
				break;
			case OPT_LOOSE:  
				os << "# opt(loose" ; 
				if( frozen.size() ) {
					os << ",modredundant";
				}
				os << ")" << endl;
				break;
			case OPT_TIGHT:  
				os << "# opt(tight" ; 
				if( frozen.size() ) {
					os << ",modredundant";
				}
				os << ")" << endl;
				break;
			default:
				return false;
		}
		

		os << "# symmetry=None" << endl ;


		if( espgrid ) {
			os << "# prop=(read,field)"  << endl;
		}


		os << endl;
		os << "MOL" << endl;
		os << endl;

		os << qi << " " << mult << endl;

		FOR_ATOMS_OF_MOL(a, mol) {
			double *c = a->GetCoordinate();
			sprintf( line, "%s %f %f %f", etab.GetSymbol(a->GetAtomicNum()), c[0], c[1], c[2] ); 
			os <<  line << endl;
		}


		os << endl;

		if( frozen.size() ) {
			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {
				if( i->size() == 2 ) {
					os << "B " << (*i)[0] << " " << (*i)[1] << " F" << endl ;
				}
			}

			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {
				if( i->size() == 3 ) {
					os << "A " << (*i)[0] << " " << (*i)[1] << " " <<  (*i)[2] <<" F" << endl ;
				}
			}

			for( std::vector< std::vector<int> >::iterator i = frozen.begin(); i!=frozen.end(); i++ ) {

				if( i->size() == 4 ) {
					os << "D " << (*i)[0] << " " << (*i)[1] << " " <<  (*i)[2] << " " << (*i)[3] << " F" << endl ;
				}
			}

		}
		else {
			os << endl;
		}

		if( espgrid ) {
			os << "@grid.dat /N" << endl;
			os << endl;
		}

		os.flush();

		return true;
	}
}

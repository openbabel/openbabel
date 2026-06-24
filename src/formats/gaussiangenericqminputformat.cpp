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

				OpenBabel::OBConversion::RegisterFormat( "gjf", this );
				OpenBabel::OBConversion::RegisterFormat( "gjc", this );
				OpenBabel::OBConversion::RegisterFormat( "gau", this );
				OpenBabel::OBConversion::RegisterFormat( "com", this, "chemical/x-gaussian-input" );

				OBConversion::RegisterOptionParam("b" , this, 0, OBConversion::OUTOPTIONS);
				OBConversion::RegisterOptionParam("u" , this, 0, OBConversion::OUTOPTIONS);
				OBConversion::RegisterOptionParam("k", this, 1, OBConversion::OUTOPTIONS);
				OBConversion::RegisterOptionParam("f", this, 1, OBConversion::OUTOPTIONS);    


			}

			virtual const char* Description() {
				string ss = string( GenericQMInputFormat::Description() );
				stringstream o;
				o << ss;
				o << "  k <string>           Keywords to include (overrides T,B,E,O,F)\n";
				o << "  f <file>             Keywords to include (can specify k and f)\n";
				o << "  b                    Include bonding information\n";
				o << "  u                    Include crystallographic unit cell, if present\n\n";
				return strdup( o.str().c_str() );
			}

			virtual const char* SpecificationURL() {

				return "http://www.gaussian.com/g_tech/g09ur.htm";
			}

			virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );
	};


	// Global variable used to register GaussianGQMInput format.
	GaussianGQMInputFormat theGaussianGenericQMInputFormat;

	static void write_unit_cell( OBMol &mol, ostream &ofs ) {

		char buffer[BUFF_SIZE];
		// Translation vectors
		OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
		if (uc ) {
			uc->FillUnitCell(&mol); // complete the unit cell with symmetry-derived atoms

			vector<vector3> cellVectors = uc->GetCellVectors();
			for (vector<vector3>::iterator i = cellVectors.begin(); i != cellVectors.end(); ++i) {
				snprintf(buffer, BUFF_SIZE, "TV       %10.5f      %10.5f      %10.5f",
						i->x(),
						i->y(),
						i->z());
				ofs << buffer << '\n';
			}
		}

	}

	static void write_bonds( OBMol &mol, ostream &ofs ) {
		char buffer[BUFF_SIZE];

		// Bonds, contributed by Daniel Mansfield
		// first, make begin.GetIdx < end.GetIdx
		OBBond* bond;
		OBAtom *atom;
		vector<OBEdgeBase*>::iterator j;
		vector<OBNodeBase*>::iterator i;
		OBAtom *bgn, *end;
		for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
		{
			if (bond->GetBeginAtomIdx() > bond->GetEndAtomIdx()) {
				bgn = bond->GetBeginAtom();
				end = bond->GetEndAtom();
				bond->SetBegin(end);
				bond->SetEnd(bgn);
			}
		}

		// this seems inefficient -- perhaps using atom neighbor iterators?
		// -GRH
		for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
		{
			ofs << endl << atom->GetIdx() << " ";
			for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
			{
				if (bond->GetBeginAtomIdx() == atom->GetIdx()) {
					snprintf(buffer, BUFF_SIZE, "%d %1.1f ", bond->GetEndAtomIdx(), (float) bond->GetBondOrder());
					ofs << buffer;
				}
			}
		} // iterate through atoms
	}

	//------------------------------------------------------------------------------
	bool GaussianGQMInputFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
	{
		char line[80];
		if( ! ParseOptions( pOb, pConv ) ) {
			return false;
		}

		bool writeUnitCell = (NULL != pConv->IsOption("u", OBConversion::OUTOPTIONS));
		bool has_keywordss  = (pConv->IsOption("k",OBConversion::OUTOPTIONS) != NULL);
		bool has_keywordsf  = (pConv->IsOption("f",OBConversion::OUTOPTIONS) != NULL);
		const char *keywordss = pConv->IsOption("k",OBConversion::OUTOPTIONS);
		const char *keywordsf = pConv->IsOption("f",OBConversion::OUTOPTIONS);


		OBMol* mol = dynamic_cast< OBMol* >(pOb);
		if( mol == 0 ) return false;

		ostream& os = *pConv->GetOutStream();

		double q = 0.;

		FOR_ATOMS_OF_MOL(atom, mol) {
			q += atom->GetPartialCharge();
		}

		int qi = (int) round(q);

		if( charge_set ) { qi = charge; }
		if( !mult_set )  { mult = mol->GetTotalSpinMultiplicity(); }


		os << std::setprecision(10);
		os << "%nprocshared="  << ncpus << endl;
		os << "%mem="  << mem << "GB" << endl;
		os << "%chk = chk" << endl;
		os << "%nosave" << endl;
		os << endl;

		// the route

		if( ! has_keywordss && ! has_keywordsf ) {

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
		}
		else {
			// If keywords are speficied, use those instead of 
			// the route assembled from parameters
			if ( has_keywordss ) {
				os << keywordss << endl;
			}
			if ( has_keywordsf ) {
				ifstream kfstream(keywordsf);
				string keyBuffer;
				if (kfstream)
				{
					while (getline(kfstream, keyBuffer))
						os << keyBuffer << endl;
				}
			}
		}




		os << endl;
		os << "MOL" << endl;
		os << endl;

		os << qi << " " << mult << endl;

		FOR_ATOMS_OF_MOL(a, mol) {
			double *c = a->GetCoordinate();
			if( a->GetIsotope() == 0 ) {
				sprintf( line, "%s %f %f %f", etab.GetSymbol(a->GetAtomicNum()), c[0], c[1], c[2] ); 
			}
			else {
				sprintf( line, "%s(Iso=%d) %f %f %f", etab.GetSymbol(a->GetAtomicNum()), a->GetIsotope(), c[0], c[1], c[2] ); 
			}
			os <<  line << endl;
		}
		if( writeUnitCell ) {
			write_unit_cell( *mol, os );
		}
		if( pConv->IsOption( "b", OBConversion::OUTOPTIONS)) {
			write_bonds( *mol, os );
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

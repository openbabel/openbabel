/**********************************************************************
	fromfile.cpp - A OBChargeModel to apply charges specified in a file

	Copyright (C) 2015 M J Harvey, Acellera Ltd
	m.j.harvey (at) acellera.com

	This file is part of the Open Babel project.
	For more information, see <http://openbabel.org/>

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
 ***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/oberror.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <openbabel/molchrg.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>


using namespace std;

namespace OpenBabel
{

	class FromFileCharges : public OBChargeModel
	{
		public:
			FromFileCharges(const char* ID) : OBChargeModel(ID, false){
			};

			const char* Description(){ return "Assign charges from file containing {'atom-name', charge} pairs"; }

			bool ComputeCharges(OBMol &mol, const char *arg );

	};

	/////////////////////////////////////////////////////////////////
	FromFileCharges theFromFileCharges("fromfile"); //Global instance

	bool read_file( const char *file, std::map< std::string, double> &q_by_name ) {
		char name[17];
		double q;

		FILE *fin = fopen( file, "r" );

		if( !fin ) {
			stringstream msg;
			msg << "Cannot open file " << file << endl;
			obErrorLog.ThrowError(__FUNCTION__, msg.str(), obError);
			return false;
		}
		while( 2 == fscanf( fin, "%16s %lf\n", name, &q ) ) {
			q_by_name.insert( std::pair<std::string, double>( string(name), q) );
		}
		fclose( fin );

		return true;
	}

	/////////////////////////////////////////////////////////////////

	bool FromFileCharges::ComputeCharges(OBMol &mol, const char *arg )
	{

		if( !arg ) {
			stringstream msg;
			msg << "Charge file argument required:" << endl
				<< "\tbabel --partialcharge fromfile:/path/to/file"<< endl
				<< "File format is one 'atom-name charge' pair per line, eg:"<<endl
				<< "\tC1\t1.0" << endl
				<< "\tO2\t-1.5" << endl;
			obErrorLog.ThrowError(__FUNCTION__, msg.str(), obError);
			return false;
		}

		std::map< std::string,  double> q_by_name;
		if( !read_file( arg, q_by_name ) ) {
			return false;
		}


		mol.SetPartialChargesPerceived();


		for ( int i = 1; i <= mol.NumAtoms(); i++)
		{
			OBAtom *a = mol.GetAtom(i);

			OBResidue *res;
			double q   = 0.;
			bool found = false;
			char *name = NULL;

			// First try atom type name
			if ( (res = a->GetResidue()) != 0 )
			{
				char *f = name  = (char*)res->GetAtomID( a ).c_str();
				for( int j = strlen(f)-1; j>=0; j-- ) { if( f[j]==' ' ){ f[j]='\0'; } } // trim trailing whitespace
				std::string ff = string(f);
				if( q_by_name.count( ff ) ) {
					q = q_by_name[ string(ff) ];
					found = true;
				}
			}
			// Then try the element symbol
			if( !found ) {
				std::string ff  = string( OBElements::GetSymbol(a->GetAtomicNum()) );
				if( q_by_name.count( ff ) ) {
					q = q_by_name[ string(ff) ];
					found = true;
				}
			}
			// Finally "*" wildcard
			if( !found ) {
				std::string ff  = string("*");
				if( q_by_name.count( "*" ) ) {
					q = q_by_name[ string(ff) ];
					found = true;
				}
			}

			if( !found ) {
				stringstream msg;
				msg << "Charge mapping for atom # " << i ;
				if( name ) {
					msg << " (" << name <<") ";
				}
				msg << "not found " <<  endl;  
				obErrorLog.ThrowError(__FUNCTION__, msg.str(), obError);
				return false;
			}

			a->SetPartialCharge( q );

		}

		OBPairData *dp = new OBPairData;
		dp->SetAttribute("PartialCharges");
		dp->SetValue("User Charges");
		dp->SetOrigin(perceived);
		mol.SetData(dp);

		OBChargeModel::FillChargeVectors(mol);

		return true;
	}

}//namespace

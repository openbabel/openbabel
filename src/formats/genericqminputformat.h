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


#ifndef __GENERICQMINPUT_H 
#define __GENERICQMINPUT_H 1

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <math.h>

using namespace std;

namespace OpenBabel
{


	//==============================================================================
	/// Class to output a point cloud on a surface around a molecule

	enum OPT_t  { OPT_NONE, OPT_LOOSE, OPT_NORMAL, OPT_TIGHT  };

	class GenericQMInputFormat : public OpenBabel::OBMoleculeFormat
	{
		protected:


			std::map<string, string> basis_sets;
			std::map<string, string> theory_sets;
			string basis_default;
			string theory_default;
			const char *format_type;


			// the parsed arguments are stored in these
			int  mult     ; //= 1;
			bool espgrid  ; //= false;
			int  mem      ; //= 2;
			int  ncpus    ; //= 1;
			std::vector<std::vector<int> > frozen;
			string  basis ; // =  basis_default;
			string  theory;
			int  opt      ; //= OPT_NONE;
			int  charge   ;  
			bool charge_set;
			bool mult_set;


		public:
			GenericQMInputFormat( const char *format_name, const char *type );
			/// Return description.
			virtual const char* Description() ;

			/// Return read/write flag: read only.
			virtual unsigned int Flags();

			/// Skip to object: used for multi-object file formats.
			virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv );

			/// Read: always return false.
			virtual bool ReadMolecule( OpenBabel::OBBase*, OpenBabel::OBConversion* );

			/// Write.
			//	virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* );

			bool ParseOptions( OBBase *pOb, OBConversion *pConv ) ;
	};

}
#endif

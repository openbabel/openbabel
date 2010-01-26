/**********************************************************************
Copyright (C) 2009 by Jeremy W. Murphy

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


/*
	Some comments about MNA further to the article, clarified via email
	with Vladimir Poroikov.

	1)	Chain (acylic) atoms come last in the sorting order.

	2)	The chain marker is exactly what it says it is.


	Some minor things that I was uncertain about in the design:

		MNAcmp() - sort() could not access it when it was a member of
		class MNAFormat... should I have friended sort() in MNAFormat?
*/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>

using namespace std;


typedef enum
{
	MNAZ_ERROR, MNAZ_H, MNAZ_C, MNAZ_N, MNAZ_O, MNAZ_F, MNAZ_Si, MNAZ_P, MNAZ_S, MNAZ_Cl, MNAZ_Ca, MNAZ_As, MNAZ_Se, MNAZ_Br, MNAZ_Li, MNAZ_B, MNAZ_Mg, MNAZ_Sn, MNAZ_Te, MNAZ_I, MNAZ_Os, MNAZ_Sc, MNAZ_Fe, MNAZ_Co, MNAZ_Sr, MNAZ_Pd, MNAZ_Be, MNAZ_K, MNAZ_V, MNAZ_Ni, MNAZ_In, MNAZ_Al, MNAZ_R
} MNAZ;


namespace OpenBabel
{

	class MNAFormat : public OBMoleculeFormat
				// Derive directly from OBFormat for objects which are not molecules.
	{
		public:
			//Register this format type ID in the constructor
			MNAFormat()
			{
				OBConversion::RegisterFormat("mna", this);
				OBConversion::RegisterOptionParam(levels_option, this, 1);
			}

			virtual const char* Description() //required
			{
				stringstream ss;

				ss << "Multilevel Neighbourhoods of Atoms (MNA)\n"
				"Iteratively generated 2D descriptors suitable for QSAR\n"
				"\n"
				"Write Options e.g. -x" << levels_option << "1 \n"
				"  " << levels_option << "# MNA Levels default " << levels << "\n";

				static const string s(ss.str());

				return s.c_str();
			};

			//Optional URL where the file format is specified
			virtual const char* SpecificationURL()
			{
				return "http://pubs.acs.org/doi/abs/10.1021/ci980335o\n"
               "http://ebook.rsc.org/?DOI=10.1039/9781847558879";
			};


			virtual unsigned int Flags()
			{
				return NOTREADABLE;	// possibly WRITEONEONLY??
			};


			////////////////////////////////////////////////////
			/// Declarations for the "API" interface functions. Definitions are below
			// virtual bool ReadMolecule ( OBBase* pOb, OBConversion* pConv );
			virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
			static MNAZ Z2MNAZ(unsigned);


		private:
			string MNAize(OBAtom * const, unsigned);

			static string const chain_marker;
			static string const open_nbor_list;
			static string const close_nbor_list;
			static const char * const levels_option;
			static unsigned levels;
	};


	string const MNAFormat::chain_marker = "-";
	string const MNAFormat::open_nbor_list = "(";
	string const MNAFormat::close_nbor_list = ")";
	const char * const MNAFormat::levels_option = "L";
	unsigned MNAFormat::levels = 2; //default
	////////////////////////////////////////////////////

	//Make an instance of the format class
	MNAFormat theMNAFormat;

	////////////////////////////////////////////////////////////////

	bool MNAcmp(OBAtom const * const, OBAtom const * const);


	bool MNAFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
	{
		OBMol* pmol = dynamic_cast<OBMol*>(pOb);

		if (pmol == NULL)
			return false;

		ostream& ofs = *pConv->GetOutStream();

		/** Write the representation of the OBMol molecule to the output stream **/

		// To find out whether this is the first molecule to be output...

		if (pConv->GetOutputIndex() == 1)
			ofs << "# The contents of this file were derived from " << pConv->GetInFilename() << endl;

		//To use an output option
		if (!pConv->IsOption("n"))      //OBConversion::OUTOPTIONS is the default
			ofs << "# Title = " << pmol->GetTitle() << endl;

		const char * const p = pConv->IsOption(levels_option, OBConversion::OUTOPTIONS);

		if (p)
		{
			stringstream ss(p);

			ss >> levels;	// 0 on error

			if (levels > 10000)
			{
				cerr << "Levels = " << levels << " will almost certainly crash and was probably a mistake." << endl;
				return false;
			}
			else
				if (levels > 10)
					cerr << "Levels > 10 is probably not very useful.  (If it is, let me know!)" << endl;
		}

		OBMol &mol = *pmol;

		vector<OBAtom*>::iterator i;
		OBAtom *atom;
		mol.AddHydrogens();

		for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
			ofs << MNAize(atom, levels) << endl;

		ofs << endl;

		// ... or the last
		if (!pConv->IsLast())
			ofs << "$$$$" << endl;

		return true; //or false to stop converting
	}


	string MNAFormat::MNAize(OBAtom * const atom, unsigned level)
	{
		string s;

		// Append level 0 string for this item, e.g.: "-H"
		s += (atom->IsInRing() ? "" : chain_marker) + etab.GetSymbol(atom->GetAtomicNum());

		if (level != 0)
		{
			// Possible problem: what if it's a single-atom molecule?  Unlikely, but no-one likes a segfault.
			// Get the list of neighbours and re-order them.
			vector<OBAtom *> nbors;
			FOR_NBORS_OF_ATOM(nbor, atom)
				nbors.push_back(&*nbor);
			sort(nbors.begin(), nbors.end(), MNAcmp);

			s += open_nbor_list;

			// Recurse into neighbours

			for (vector<OBAtom *>::iterator nbor = nbors.begin(); nbor != nbors.end(); ++nbor)
				s += MNAize(*nbor, level - 1);

			s += close_nbor_list;
		}

		return s;
	}


	// Cyclic atoms come before acyclic, then in MNA type order
	bool MNAcmp(OBAtom const * const A, OBAtom const * const B)
	{
		bool altb;


		if (A->IsInRing())
		{
			if (B->IsInRing())
			{
				if (MNAFormat::Z2MNAZ(A->GetAtomicNum()) < MNAFormat::Z2MNAZ(B->GetAtomicNum()))
					altb = true;
				else
					altb = false; // Will this need to recurse here?  (Or does atomic numbering take care of that?)
			}
			else
				altb = true;
		}
		else
		{
			if (B->IsInRing())
				altb = false;
			else
			{
				if (MNAFormat::Z2MNAZ(A->GetAtomicNum()) < MNAFormat::Z2MNAZ(B->GetAtomicNum()))
					altb = true;
				else
					altb = false; // And here?
			}
		}

		return altb;
	}


	// Big, boring look-up table to establish the MNA type.
	// It should also be noted that in their article, Filimonov et al refer to element Jl, which is Jolotium, a former name for Dubnium.
	MNAZ MNAFormat::Z2MNAZ(unsigned int Z)
	{
		MNAZ type;

		switch (Z)
		{
			case 0: type = MNAZ_ERROR; break;
			case 1: type = MNAZ_H; break;
			case 2: type = MNAZ_R; break;
			case 3: type = MNAZ_Li; break;
			case 4: type = MNAZ_Be; break;
			case 5: type = MNAZ_B; break;
			case 6: type = MNAZ_C; break;
			case 7: type = MNAZ_N; break;
			case 8: type = MNAZ_O; break;
			case 9: type = MNAZ_F; break;
			case 10: type = MNAZ_R; break;
			case 11: type = MNAZ_Li; break;
			case 12: type = MNAZ_Mg; break;
			case 13: type = MNAZ_Al; break;
			case 14: type = MNAZ_Si; break;
			case 15: type = MNAZ_P; break;
			case 16: type = MNAZ_S; break;
			case 17: type = MNAZ_Cl; break;
			case 18: type = MNAZ_R; break;
			case 19: type = MNAZ_K; break;
			case 20: type = MNAZ_Ca; break;
			case 21: type = MNAZ_Sc; break;
			case 22: type = MNAZ_Sc; break;
			case 23: type = MNAZ_V; break;
			case 24: type = MNAZ_V; break;
			case 25: type = MNAZ_Mg; break;
			case 26: type = MNAZ_Fe; break;
			case 27: type = MNAZ_Co; break;
			case 28: type = MNAZ_Ni; break;
			case 29: type = MNAZ_Ni; break;
			case 30: type = MNAZ_Be; break;
			case 31: type = MNAZ_Al; break;
			case 32: type = MNAZ_Ni; break;
			case 33: type = MNAZ_As; break;
			case 34: type = MNAZ_Se; break;
			case 35: type = MNAZ_Br; break;
			case 36: type = MNAZ_R; break;
			case 37: type = MNAZ_K; break;
			case 38: type = MNAZ_Sr; break;
			case 39: type = MNAZ_Al; break;
			case 40: type = MNAZ_Sc; break;
			case 41: type = MNAZ_V; break;
			case 42: type = MNAZ_V; break;
			case 43: type = MNAZ_V; break;
			case 44: type = MNAZ_Ni; break;
			case 45: type = MNAZ_Ni; break;
			case 46: type = MNAZ_Pd; break;
			case 47: type = MNAZ_Ni; break;
			case 48: type = MNAZ_Be; break;
			case 49: type = MNAZ_In; break;
			case 50: type = MNAZ_Sn; break;
			case 51: type = MNAZ_Co; break;
			case 52: type = MNAZ_Te; break;
			case 53: type = MNAZ_I; break;
			case 54: type = MNAZ_R; break;
			case 55: type = MNAZ_K; break;
			case 56: type = MNAZ_Sr; break;
			case 57: type = MNAZ_In; break;
			case 58: type = MNAZ_In; break;
			case 59: type = MNAZ_In; break;
			case 60: type = MNAZ_In; break;
			case 61: type = MNAZ_In; break;
			case 62: type = MNAZ_In; break;
			case 63: type = MNAZ_In; break;
			case 64: type = MNAZ_Al; break;
			case 65: type = MNAZ_Al; break;
			case 66: type = MNAZ_Al; break;
			case 67: type = MNAZ_Al; break;
			case 68: type = MNAZ_Al; break;
			case 69: type = MNAZ_Al; break;
			case 70: type = MNAZ_Al; break;
			case 71: type = MNAZ_Al; break;
			case 72: type = MNAZ_Fe; break;
			case 73: type = MNAZ_Fe; break;
			case 74: type = MNAZ_Co; break;
			case 75: type = MNAZ_B; break;
			case 76: type = MNAZ_Os; break;
			case 77: type = MNAZ_Os; break;
			case 78: type = MNAZ_Pd; break;
			case 79: type = MNAZ_Pd; break;
			case 80: type = MNAZ_Be; break;
			case 81: type = MNAZ_Al; break;
			case 82: type = MNAZ_Sn; break;
			case 83: type = MNAZ_Ni; break;
			case 84: type = MNAZ_Te; break;
			case 85: type = MNAZ_I; break;
			case 86: type = MNAZ_R; break;
			case 87: type = MNAZ_K; break;
			case 88: type = MNAZ_Sr; break;
			case 89: type = MNAZ_R; break;
			case 90: type = MNAZ_R; break;
			case 91: type = MNAZ_R; break;
			case 92: type = MNAZ_R; break;
			case 93: type = MNAZ_R; break;
			case 94: type = MNAZ_R; break;
			case 95: type = MNAZ_R; break;
			case 96: type = MNAZ_R; break;
			case 97: type = MNAZ_R; break;
			case 98: type = MNAZ_R; break;
			case 99: type = MNAZ_R; break;
			case 100: type = MNAZ_R; break;
			case 101: type = MNAZ_R; break;
			case 102: type = MNAZ_R; break;
			case 103: type = MNAZ_R; break;
			case 104: type = MNAZ_ERROR; break;
			case 105: type = MNAZ_R; break;
			case 106: type = MNAZ_R; break;
			default: type = MNAZ_ERROR; break;
		}

		return type;
	}
}

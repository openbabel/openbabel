/**********************************************************************
Copyright (C) 2010 by Jeremy W. Murphy

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 or later of the License.

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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/ring.h>

using namespace std;


typedef enum
{
	MNAZ_ERROR, MNAZ_H, MNAZ_C, MNAZ_N, MNAZ_O, MNAZ_F, MNAZ_Si, MNAZ_P, MNAZ_S, MNAZ_Cl, MNAZ_Ca, MNAZ_As, MNAZ_Se, MNAZ_Br, MNAZ_Li, MNAZ_B, MNAZ_Mg, MNAZ_Sn, MNAZ_Te, MNAZ_I, MNAZ_Os, MNAZ_Sc, MNAZ_Fe, MNAZ_Co, MNAZ_Sr, MNAZ_Pd, MNAZ_Be, MNAZ_K, MNAZ_V, MNAZ_Ni, MNAZ_In, MNAZ_Al
} MNAZ;


namespace OpenBabel
{

	class MNAFormat : public OBMoleculeFormat
	{
		public:
			MNAFormat()
			{
				OBConversion::RegisterFormat("mna", this);
				OBConversion::RegisterOptionParam(levels_option, this, 1);
			}

			virtual const char* Description()
			{
        stringstream ss;
				ss <<

"Multilevel Neighborhoods of Atoms (MNA)\n"
"Iteratively generated 2D descriptors suitable for QSAR\n"
"Multilevel Neighborhoods of Atoms (MNA) descriptors are\n"
"2D molecular fragments suitable for use in QSAR modelling [fpbg99]_.\n"
"The format outputs a complete descriptor fingerprint per\n"
"molecule. Thus, a 27-atom (including hydrogen) molecule would\n"
"result in 27 descriptors, one per line.\n\n"

"MNA descriptors are generated recursively. Starting at the origin,\n"
"each atom is appended to the descriptor immediately followed by a\n"
"parenthesized list of its neighbours. This process iterates until the\n"
"specified distance from the origin, also known as the depth of the\n"
"descriptor.\n\n"

"Elements are simplified into 32 groups. Each group has a representative\n"
"symbol used to stand for any element in that group:\n\n"

"==== ========\n"
"Type Elements\n"
"==== ========\n"
"H    H\n"
"C    C\n"
"N    N\n"
"O    O\n"
"F    F\n"
"Si   Si\n"
"P    P\n"
"S    S\n"
"Cl   Cl\n"
"Ca   Ca\n"
"As   As\n"
"Se   Se\n"
"Br   Br\n"
"Li   Li, Na\n"
"B    B, Re\n"
"Mg   Mg, Mn\n"
"Sn   Sn, Pb\n"
"Te   Te, Po\n"
"I    I, At\n"
"Os   Os, Ir\n"
"Sc   Sc, Ti, Zr\n"
"Fe   Fe, Hf, Ta\n"
"Co   Co, Sb, W\n"
"Sr   Sr, Ba, Ra\n"
"Pd   Pd, Pt, Au\n"
"Be   Be, Zn, Cd, Hg\n"
"K    K, Rb, Cs, Fr\n"
"V    V, Cr, Nb, Mo, Tc\n"
"Ni   Ni, Cu, Ge, Ru, Rh, Ag, Bi\n"
"In   In, La, Ce, Pr, Nd, Pm, Sm, Eu\n"
"Al   Al, Ga, Y, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Tl\n"
"R    R, He, Ne, Ar, Kr, Xe, Rn, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Db, Jl\n"
"==== ========\n\n"

"Acyclic atoms are preceded by a hyphen \"-\" mark.\n\n"

"Here's the multi-level neighborhood for the molecule\n"
"represented by the SMILES string CC(=O)Cl::\n\n"

" # The contents of this file were derived from\n"
" # Title = acid chloride\n"
" -C(-H(-C)-H(-C)-H(-C)-C(-C-O-Cl))\n"
" -C(-C(-H-H-H-C)-O(-C)-Cl(-C))\n"
" -O(-C(-C-O-Cl))\n"
" -Cl(-C(-C-O-Cl))\n"
" -H(-C(-H-H-H-C))\n"
" -H(-C(-H-H-H-C))\n"
" -H(-C(-H-H-H-C))\n\n"

".. [fpbg99] Dmitrii Filimonov, Vladimir Poroikov, Yulia Borodina, and\n"
"            Tatyana Gloriozova. **Chemical Similarity Assessment through\n"
"            Multilevel Neighborhoods of Atoms: Definition and Comparison with\n"
"            the Other Descriptors.** *J. Chem. Inf. Comput. Sci.* **1999**, *39*, 666-670.\n"
"            [`Link <https://doi.org/10.1021/ci980335o>`_]\n\n"

        "Write Options e.g. -x" << levels_option << "1 \n"
				"  " << levels_option << "#  Levels (default = " << levels << ")\n\n"
        ;

				static const string s(ss.str());

				return s.c_str();
			};

			virtual const char* SpecificationURL()
			{
				return "http://openbabel.org/wiki/Multilevel_Neighborhoods_of_Atoms";
			};

			virtual unsigned int Flags()
			{
				return NOTREADABLE;	// possibly WRITEONEONLY??
			};

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
	unsigned MNAFormat::levels = 2; // default

	MNAFormat theMNAFormat;

	////////////////////////////////////////////////////////////////

	bool MNAcmp(OBAtom const * const, OBAtom const * const);


	bool MNAFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
	{
		OBMol* pmol = dynamic_cast<OBMol*>(pOb);

		if (pmol == NULL)
			return false;

		ostream& ofs = *pConv->GetOutStream();

		if (pConv->GetOutputIndex() == 1)
			ofs << "# The contents of this file were derived from " << pConv->GetInFilename() << endl;

		if (!pConv->IsOption("n"))
			ofs << "# Title = " << pmol->GetTitle() << endl;

		const char * const p = pConv->IsOption(levels_option);

		if (p)
		{
			stringstream ss(p), error_msg;


			ss >> levels;

			if (!(ss.rdstate() & (stringstream::failbit | stringstream::badbit)))
			{


				if (levels > 10000)
				{
					error_msg << "Levels = " << levels << " will almost certainly crash and was probably a mistake." << endl;
					obErrorLog.ThrowError(__FUNCTION__, error_msg.str(), obError);
					return false;
				}
				else
					if (levels > 10)
					{
						error_msg << "Levels > 10 is probably not very useful.  (If it is, let me know!)" << endl;
						obErrorLog.ThrowError(__FUNCTION__, error_msg.str(), obWarning);
					}
			}
			else
			{
				error_msg << "Error reading levels value: " << ss.str() << endl;
				obErrorLog.ThrowError(__FUNCTION__, error_msg.str(), obError);
				return false;
			}
		}

		OBMol &mol = *pmol;

		if (pConv->IsOption("d", OBConversion::GENOPTIONS))
		{
			if (pConv->GetOutputIndex() == 1)
				obErrorLog.ThrowError(__FUNCTION__, "MNA includes hydrogens by definition, just be aware of that.", obInfo);
			ofs << "# Hydrogens deleted explicitly." << endl;
			mol.DeleteHydrogens();
		}
		else
			mol.AddHydrogens();

		FOR_ATOMS_OF_MOL(atom, mol)
			ofs << MNAize(&*atom, levels) << endl;

		if (!pConv->IsLast())
			ofs << "$$$$" << endl;

		return true;
	}


	string MNAFormat::MNAize(OBAtom * const atom, unsigned level)
	{
		string s;

		// Append level 0 string for this item, e.g.: "-H"
		s += (atom->IsInRing() ? "" : chain_marker) + OBElements::GetSymbol(atom->GetAtomicNum());

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
	MNAZ MNAFormat::Z2MNAZ(unsigned int Z)
	{
		MNAZ type;

		switch (Z)
		{
			case 0: type = MNAZ_ERROR; break;
			case 1: type = MNAZ_H; break;
			case 2: type = MNAZ_ERROR; break;
			case 3: type = MNAZ_Li; break;
			case 4: type = MNAZ_Be; break;
			case 5: type = MNAZ_B; break;
			case 6: type = MNAZ_C; break;
			case 7: type = MNAZ_N; break;
			case 8: type = MNAZ_O; break;
			case 9: type = MNAZ_F; break;
			case 10: type = MNAZ_ERROR; break;
			case 11: type = MNAZ_Li; break;
			case 12: type = MNAZ_Mg; break;
			case 13: type = MNAZ_Al; break;
			case 14: type = MNAZ_Si; break;
			case 15: type = MNAZ_P; break;
			case 16: type = MNAZ_S; break;
			case 17: type = MNAZ_Cl; break;
			case 18: type = MNAZ_ERROR; break;
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
			case 36: type = MNAZ_ERROR; break;
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
			case 54: type = MNAZ_ERROR; break;
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
			case 86: type = MNAZ_ERROR; break;
			case 87: type = MNAZ_K; break;
			case 88: type = MNAZ_Sr; break;
			case 89: type = MNAZ_ERROR; break;
			case 90: type = MNAZ_ERROR; break;
			case 91: type = MNAZ_ERROR; break;
			case 92: type = MNAZ_ERROR; break;
			case 93: type = MNAZ_ERROR; break;
			case 94: type = MNAZ_ERROR; break;
			case 95: type = MNAZ_ERROR; break;
			case 96: type = MNAZ_ERROR; break;
			case 97: type = MNAZ_ERROR; break;
			case 98: type = MNAZ_ERROR; break;
			case 99: type = MNAZ_ERROR; break;
			case 100: type = MNAZ_ERROR; break;
			case 101: type = MNAZ_ERROR; break;
			case 102: type = MNAZ_ERROR; break;
			case 103: type = MNAZ_ERROR; break;
			case 104: type = MNAZ_ERROR; break;
			case 105: type = MNAZ_ERROR; break;
			case 106: type = MNAZ_ERROR; break;
			default: type = MNAZ_ERROR; break;
		}

		return type;
	}
}

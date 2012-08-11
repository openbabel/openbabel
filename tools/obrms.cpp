/**********************************************************************
 obrms = Returns the rms between two chemically identical structures
 Derived from obfit.
 *
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

/*
 Require a fixed molecule and a set of molecules to compare against.
 example of command line:
 obrms ref.sdf test.sdf
 */

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/query.h>
#include <openbabel/isomorphism.h>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <cassert>
#include <sstream>

using namespace std;
using namespace OpenBabel;

class AtomDistanceSorter
{
	vector3 ref;
public:
	AtomDistanceSorter(OBAtom *r) :
			ref(r->GetVector())
	{

	}
	bool operator()(OBAtom *l, OBAtom *r) const
	{
		double ld = ref.distSq(l->GetVector());
		double rd = ref.distSq(r->GetVector());

		return ld < rd;
	}
};
/* class to perform graph matching between two molecules.
 * Is initialized with the reference molecule.
 * Will figure out the atom correspondences and compute the rmsd between the
 * ref mol and a provided test mol.
 */
class Matcher
{
	const OBMol& ref;
	OBQuery *query;
	OBIsomorphismMapper *mapper;

	class MapRMSDFunctor: public OBIsomorphismMapper::Functor
	{
	private:
		const OBMol& ref;
		const OBMol& test;
		double bestRMSD;
	public:
		MapRMSDFunctor(const OBMol& r, const OBMol& t) :
				ref(r), test(t), bestRMSD(HUGE_VAL)
		{
		}

		bool operator()(OBIsomorphismMapper::Mapping &map)
		{
			double rmsd = 0;
			for (unsigned i = 0, n = map.size(); i < n; i++)
			{
				//obmol indices are 1-indexed while the mapper is zero indexed
				const OBAtom *ratom = ref.GetAtom(map[i].first + 1);
				const OBAtom *tatom = test.GetAtom(map[i].second + 1);
				assert(ratom && tatom);
				vector3 rvec = ratom->GetVector();
				vector3 tvec = tatom->GetVector();
				rmsd += rvec.distSq(tvec);
			}
			rmsd /= map.size();
			rmsd = sqrt(rmsd);
			if (rmsd < bestRMSD)
				bestRMSD = rmsd;
			// check all possible mappings
			return false;
		}

		double getRMSD() const
		{
			return bestRMSD;
		}
	};

public:
	Matcher(OBMol& mol) :
			ref(mol), query(NULL), mapper(NULL)
	{
		query = CompileMoleculeQuery(&mol);
		mapper = OBIsomorphismMapper::GetInstance(query);
	}

	~Matcher()
	{
		if (query)
			delete query;
		if (mapper)
			delete mapper;
	}

	//computes a correspondence between the ref mol and test (exhaustively)
	//and returns the rmsd; returns infinity if unmatchable
	double computeRMSD(OBMol& test)
	{
		MapRMSDFunctor funct(ref, test);

		mapper->MapGeneric(funct, &test);
		return funct.getRMSD();
	}
};

//preprocess molecule into a standardized state for heavy atom rmsd computation
static void processMol(OBMol& mol)
{
	//isomorphismmapper wants isomorphic atoms to have the same aromatic and ring state,
	//but these proporties aren't reliable enough to be trusted in evaluating molecules
	//should be considered the same based solely on connectivity
	
	mol.DeleteHydrogens(); //heavy atom rmsd
	for(OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); aitr++)
	{
		OBAtom *a = *aitr;
		a->UnsetAromatic();
		a->SetInRing();
	}
	for(OBBondIterator bitr = mol.BeginBonds(); bitr != mol.EndBonds(); bitr++)
	{
		OBBond *b = *bitr;
		b->UnsetAromatic();
		b->SetBondOrder(1);
		b->SetInRing();
	}
	//avoid recomputations
	mol.SetHybridizationPerceived();
	mol.SetRingAtomsAndBondsPerceived();
	mol.SetAromaticPerceived();
}
///////////////////////////////////////////////////////////////////////////////
//! \brief compute rms between chemically identical molecules
int main(int argc, char **argv)
{
	bool firstOnly = false;
	if (argc != 3 && argc != 4)
	{
		cerr << "Usage: " << argv[0]
				<< " [-firstonly] <reference structure(s)> <comparison structure(s)>\n";
		cerr << "Computes the heavy-atom RMSD of identical compound structures.\n";
		cerr << "Structures in multi-structure files are compared one-by-one unless -firstonly\n" 
		<< "is passed, in which case only the first structure in the reference file is used.\n";
		exit(-1);
	}

	char *fileRef = argv[1];
	char *fileTest = argv[2];

	if (argc == 4)
	{
		//if iterate is passed as first command, try to match structures in first file to strucutres in second
		if (strcmp("-firstonly", argv[1]) != 0)
		{
			cerr << "Usage: " << argv[0]
					<< " [-firstonly] <reference structure(s)> <comparison structure(s)>\n";
			exit(-1);
		}

		fileRef = argv[2];
		fileTest = argv[3];
		firstOnly = true;
	}

	//open mols
	OBConversion refconv;
	OBFormat *refFormat = refconv.FormatFromExt(fileRef);
	if (!refFormat || !refconv.SetInFormat(refFormat)
			|| !refconv.SetOutFormat("SMI"))
	{
		cerr << "Cannot read reference molecule format!" << endl;
		exit(-1);
	}

	OBConversion testconv;
	OBFormat *testFormat = testconv.FormatFromExt(fileTest);
	if (!testFormat || !testconv.SetInAndOutFormats(testFormat, testFormat))
	{
		cerr << "Cannot read reference molecule format!" << endl;
		exit(-1);
	}

	//read reference
	ifstream ifsref;
	OBMol molref;

	ifsref.open(fileRef);
	if (!ifsref)
	{
		cerr << "Cannot read fixed molecule file: " << fileRef << endl;
		exit(-1);
	}

	//check comparison file
	ifstream ifstest;
	ifstest.open(fileTest);
	if (!ifstest)
	{
		cerr << "Cannot read file: " << fileTest << endl;
		exit(-1);
	}

	while (refconv.Read(&molref, &ifsref))
	{
		processMol(molref);
		Matcher matcher(molref);// create the matcher
		OBMol moltest;
		while (testconv.Read(&moltest, &ifstest))
		{
			if (moltest.Empty())
				break;

			processMol(moltest);

			double rmsd = matcher.computeRMSD(moltest);

			cout << "RMSD " << moltest.GetTitle() << " " << rmsd << "\n";
			if (!firstOnly)
			{
				break;
			}
		}
	}
	return (0);
}

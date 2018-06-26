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

#undef _GLIBCXX_DEBUG /* unless you have boost libraries built with this, you do _not_ want this defined */
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

#include <sstream>
#include <memory>
#include <boost/unordered_map.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

using namespace std;
using namespace boost;
using namespace boost::program_options;
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
	std::shared_ptr<OBQuery> query;
	std::shared_ptr<OBIsomorphismMapper> mapper;

	class MapRMSDFunctor: public OBIsomorphismMapper::Functor
	{
	private:
		const OBMol& ref;
		OBMol& test;
		double bestRMSD;
		bool minimize;
		public:
		//will modify t if min is true
		MapRMSDFunctor(const OBMol& r, OBMol& t, bool min = false) :
				ref(r), test(t), bestRMSD(HUGE_VAL), minimize(min)
		{
		}

		bool operator()(OBIsomorphismMapper::Mapping &map)
		{
			unsigned N = map.size();
			double refcoord[N * 3];
			double testcoord[N * 3];

			for (unsigned i = 0; i < N; i++)
			{
				//obmol indices are 1-indexed while the mapper is zero indexed
				const OBAtom *ratom = ref.GetAtom(map[i].first + 1);
				const OBAtom *tatom = test.GetAtom(map[i].second + 1);
				assert(ratom && tatom);

				for (unsigned c = 0; c < 3; c++)
				{
					refcoord[3 * i + c] = ratom->GetVector()[c];
					testcoord[3 * i + c] = tatom->GetVector()[c];
				}
			}

			if (minimize)
			{
				double rmatrix[3][3] =
						{ 0 };

				double rave[3] =
						{ 0, 0, 0 };
				double tave[3] =
						{ 0, 0, 0 };
				//center
				for (unsigned i = 0; i < N; i++)
				{
					for (unsigned c = 0; c < 3; c++)
					{
						rave[c] += refcoord[3 * i + c];
						tave[c] += testcoord[3 * i + c];
					}
				}

				for (unsigned c = 0; c < 3; c++)
				{
					rave[c] /= N;
					tave[c] /= N;
				}

				for (unsigned i = 0; i < N; i++)
				{
					for (unsigned c = 0; c < 3; c++)
					{
						refcoord[3 * i + c] -= rave[c];
						testcoord[3 * i + c] -= tave[c];
					}
				}

				qtrfit(refcoord, testcoord, N, rmatrix);
				rotate_coords(testcoord, rmatrix, N);
				
				for (unsigned i = 0; i < N; i++)
				{
					//with minimize on, change coordinates
					OBAtom *tatom = test.GetAtom(map[i].second + 1);
					tatom->SetVector(testcoord[3*i]+rave[0], testcoord[3*i+1]+rave[1], testcoord[3*i+2]+rave[2]);
				}
			}

			double rmsd = calc_rms(refcoord, testcoord, N);

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
	Matcher(OBMol& mol) : ref(mol)
	{
		query = std::shared_ptr<OBQuery>(CompileMoleculeQuery(&mol));
		mapper = std::shared_ptr<OBIsomorphismMapper>(OBIsomorphismMapper::GetInstance(query.get()));
	}


	//computes a correspondence between the ref mol and test (exhaustively)
	//and returns the rmsd; returns infinity if unmatchable
	double computeRMSD(OBMol& test, bool minimize = false)
	{
		MapRMSDFunctor funct(ref, test, minimize);

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
	for (OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); aitr++)
	{
		OBAtom *a = *aitr;
		a->UnsetAromatic();
		a->SetInRing();
	}
	for (OBBondIterator bitr = mol.BeginBonds(); bitr != mol.EndBonds(); bitr++)
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
	bool minimize = false;
	bool separate = false;
	bool help = false;
	string fileRef;
	string fileTest;
	string fileOut;
	
	program_options::options_description desc("Allowed options");
	desc.add_options()
	("reference", value<string>(&fileRef)->required(),
			"reference structure(s) file")
	("test", value<string>(&fileTest)->required(), "test structure(s) file")
	("firstonly,f", bool_switch(&firstOnly),
			"use only the first structure in the reference file")
	("minimize,m", bool_switch(&minimize), "compute minimum RMSD")
  ("separate,s", bool_switch(&separate), "separate reference file into constituent molecules and report best RMSD")
	("out", value<string>(&fileOut), "re-oriented test structure output")
	("help", bool_switch(&help), "produce help message");

	positional_options_description pd;
	pd.add("reference", 1).add("test", 1);

	variables_map vm;
	try
	{
		store(
				command_line_parser(argc, argv).options(desc).positional(pd).run(),
				vm);
		notify(vm);
	} catch (boost::program_options::error& e)
	{
		std::cerr << "Command line parse error: " << e.what() << '\n' << desc
				<< '\n';
		exit(-1);
	}

	if (help)
	{
		cout
		<< "Computes the heavy-atom RMSD of identical compound structures.\n";
		cout << desc;
		exit(0);
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

	OBConversion outconv;
	OBFormat *outFormat = outconv.FormatFromExt(fileOut);
	if(fileOut.size() > 0)
	{
		if(!outFormat || !outconv.SetInAndOutFormats(outFormat, outFormat))
		{
			cerr << "Do not understand output format!" << endl;
			exit(-1);
		}
	}
	
	//read reference
	OBMol molref;
	std::ifstream uncompressed_inmol(fileRef.c_str());
	iostreams::filtering_stream<iostreams::input> ifsref;
	string::size_type pos = fileRef.rfind(".gz");
	if (pos != string::npos)
	{
		ifsref.push(iostreams::gzip_decompressor());
	}
	ifsref.push(uncompressed_inmol);

	if (!ifsref || !uncompressed_inmol)
	{
		cerr << "Cannot read fixed molecule file: " << fileRef << endl;
		exit(-1);
	}

	//check comparison file
	std::ifstream uncompressed_test(fileTest.c_str());
	iostreams::filtering_stream<iostreams::input> ifstest;
	pos = fileTest.rfind(".gz");
	if (pos != string::npos)
	{
		ifstest.push(iostreams::gzip_decompressor());
	}
	ifstest.push(uncompressed_test);

	if (!ifstest || !uncompressed_test)
	{
		cerr << "Cannot read file: " << fileTest << endl;
		exit(-1);
	}
	
	std::ofstream out;
	if(fileOut.size() > 0)
	{
		out.open(fileOut.c_str());
	}

	while (refconv.Read(&molref, &ifsref))
	{
	  vector<OBMol> refmols;
	  if(separate) {
	    refmols = molref.Separate();
	  } else {
	    refmols.push_back(molref);
	  }

	  vector<Matcher> matchers;
	  for(unsigned i = 0, n = refmols.size(); i < n; i++) {
	    processMol(refmols[i]);
	    Matcher matcher(refmols[i]); // create the matcher
	    matchers.push_back(matcher);
	  }

		OBMol moltest;
		while (testconv.Read(&moltest, &ifstest))
		{
			if (moltest.Empty())
				break;

			processMol(moltest);

			double bestRMSD = HUGE_VAL;
			for(unsigned i = 0, n = matchers.size(); i < n; i++) {
			  double rmsd = matchers[i].computeRMSD(moltest, minimize);
			  if(rmsd < bestRMSD) bestRMSD = rmsd;
			}

			cout << "RMSD " << molref.GetTitle() << ":" <<  moltest.GetTitle() << " " << bestRMSD << "\n";
			
			if(out)
			{
				outconv.Write(&moltest, &out);
			}
			if (!firstOnly)
			{
				break;
			}
		}
	}
	return (0);
}

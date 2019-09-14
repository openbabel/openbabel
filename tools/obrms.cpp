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
#include <cstdlib>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/query.h>
#include <openbabel/isomorphism.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/obutil.h>

#include "getopt.h"

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
	const OBMol *ref;
	obsharedptr<OBQuery> query;
	obsharedptr<OBIsomorphismMapper> mapper;

	class MapRMSDFunctor: public OBIsomorphismMapper::Functor
	{
	private:
		const OBMol *ref;
		OBMol& test;
		double bestRMSD;
		bool minimize;
		public:
		//will modify t if min is true
		MapRMSDFunctor(const OBMol* r, OBMol& t, bool min = false) :
				ref(r), test(t), bestRMSD(HUGE_VAL), minimize(min)
		{
		}

		bool operator()(OBIsomorphismMapper::Mapping &map)
		{
			unsigned N = map.size();
			double *refcoord = (double*)alloca(sizeof(double)*N * 3);
			double *testcoord = (double*)alloca(sizeof(double)*N * 3);

			for (unsigned i = 0; i < N; i++)
			{
				//obmol indices are 1-indexed while the mapper is zero indexed
				const OBAtom *ratom = ref->GetAtom(map[i].first + 1);
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
	Matcher(OBMol& mol) : ref(&mol)
	{
		query = obsharedptr<OBQuery>(CompileMoleculeQuery(&mol));
		mapper = obsharedptr<OBIsomorphismMapper>(OBIsomorphismMapper::GetInstance(query.get()));
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
		a->SetAromatic(false);
		a->SetInRing();
	}
	for (OBBondIterator bitr = mol.BeginBonds(); bitr != mol.EndBonds(); bitr++)
	{
		OBBond *b = *bitr;
		b->SetAromatic(false);
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
	bool docross = false;
	string fileRef;
	string fileTest;
	string fileOut;

	const char *helpmsg =
	 "obrms: Computes the heavy-atom RMSD of identical compound structures.\n"
	  "Usage: obrms reference_file [test_file]\n"
	  "Options:\n"
    "\t -o, --out        re-oriented test structure output\n"
	  "\t -f, --firstonly  use only the first structure in the reference file\n"
	  "\t -m, --minimize   compute minimum RMSD\n"
	  "\t -x, --cross      compute all n^2 RMSDs between molecules of reference file\n"
	  "\t -s, --separate   separate reference file into constituent molecules and report best RMSD\n"
	  "\t -h, --help       help message\n";
	struct option long_options[] = {
	    {"firstonly", no_argument, 0, 'f'},
	    {"minimize", no_argument, 0, 'm'},
	    {"cross", no_argument, 0, 'x'},
	    {"separate", no_argument, 0, 's'},
	    {"out", required_argument, 0, 'o'},
	    {"help", no_argument, 0, 'h'}
	};
	int option_index = 0;
	int c = 0;
	while ((c = getopt_long(argc, argv, "hfmxso:", long_options, &option_index) ) > 0) {
	  switch(c) {
	    case 'o':
	      fileOut = optarg;
	      break;
	    case 'f':
	      firstOnly = true;
	      break;
	    case 'm':
	      minimize = true;
	      break;
	    case 'x':
	      docross = true;
	      break;
	    case 's':
	      separate = true;
	      break;
	    case 'h':
	      cout << helpmsg;
	      exit(0);
	      break;
	    default:
	      cerr << "Unrecognized option: " << c << "\n";
	      exit(-1);
	  }
	}

	if(optind < argc) {
	  fileRef = argv[optind];
	  optind++;
	}
	if(optind < argc) {
	  fileTest = argv[optind];
	  optind++;
	}

	if(optind < argc) {
	  cerr << "Unrecognized argument: " << argv[optind];
	  exit(-1);
	}

	if(!docross && fileTest.size() == 0) {
    cerr << helpmsg;
	  cerr << "Command line parse error: test file is required but missing\n";
	  exit(-1);
	}

	//open mols
	OBConversion refconv(fileRef);

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

	std::ofstream out;
	if(fileOut.size() > 0)
	{
		out.open(fileOut.c_str());
	}


  //read reference
  OBMol molref;

	if(docross) {
	  //load in the entire reference file
    vector<OBMol> refmols;
    while (refconv.Read(&molref))
    {
       processMol(molref);
       refmols.push_back(molref);
    }

    for(unsigned i = 0, n = refmols.size() ; i < n; i++) {
      OBMol& ref = refmols[i];
      Matcher matcher(ref);
      cout << ref.GetTitle();
      for(unsigned j = 0; j < n; j++) {
        OBMol& moltest = refmols[j];
        double rmsd = matcher.computeRMSD(moltest, minimize);
        cout << ", " << rmsd;
      }
      cout << "\n";
    }

	} else {

	  //check comparison file
    while (refconv.Read(&molref))
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

      OBConversion testconv(fileTest);
      OBMol moltest;
      while (testconv.Read(&moltest))
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
        if (!firstOnly) //one test molecule will be read for each reference molecule
          break;
      }
      if(firstOnly) //done with first reference mol
        break;
    }
	}
	return (0);
}

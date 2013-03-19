/**********************************************************************
calcrmsd - Calculate the RMSD of conformers with respect to reference
           structures

Copyright (C) 2010 Noel O'Boyle

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/math/align.h>

#ifdef _MSC_VER
	typedef char TCHAR;
	#include "getopt.h"
#else
	#include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

class CalcRMSD
{
public:
  CalcRMSD(OBConversion &rconv, OBConversion &tconv, double rmsd_cutoff);
  void Run();
private:
  double rmsd_cutoff;
  OBConversion &rconv, &tconv;
};

int main(int argc,char **argv)
{
  char c;
  char *program_name = argv[0];
  char *iext;

  OBConversion rconv, tconv;
  OBFormat *prFormat = NULL;
  OBFormat *ptFormat = NULL;
  float rmsd_cutoff = 0.5;

  // Parse options
  while ((c = getopt(argc, argv, "a:b:r:")) != -1)
    {
#ifdef _WIN32
	    char optopt = c;
#endif
      switch (c)
        {
        case 'a':
          iext = optarg;

          prFormat = rconv.FindFormat(iext);

          if(prFormat==NULL)
            {
              cerr << program_name << ": cannot read reference format!" << endl;
              exit(-1);
            }
          break;
        case 'b':
          iext = optarg;

          ptFormat = tconv.FindFormat(iext);

          if(ptFormat==NULL)
            {
              cerr << program_name << ": cannot read test file format!" << endl;
              exit(-1);
            }
          break;
        case 'r': // set rmsd_cutoff

          c = sscanf(optarg, "%f", &rmsd_cutoff);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -r option" << endl;
              exit (-1);
            }
          break;
        case '?':
          if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        }
    }
  int index = optind;

  if (argc-index != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " [options] <reference> <testfile>\n";
      err += "Options:\n";
      err += "   -a <format> Reference file format\n";
      err += "   -b <format> Test file format\n";
      err += "   -r <rmsd>   RMSD cutoff (default 0.5 Angstrom)\n";
      cerr << err << ends;
      exit(-1);
    }

  // Set up the input and output files
  char* reference = argv[index++];
  char* testfile = argv[index];

  ifstream rfs;
  rfs.open(reference);
  if (!rfs) {
      cerr << program_name << ": cannot read reference file!" << endl;
      exit (-1);
  }
  ifstream tfs;
  tfs.open(testfile);
  if (!tfs) {
      cerr << program_name << ": cannot read test file!" << endl;
      exit (-1);
  }

  if (!prFormat) {
    prFormat = rconv.FormatFromExt(reference);
    if (!prFormat) {
      cerr << program_name << ": cannot read reference format!" << endl;
      return (-1);
    }
  }
  if (!ptFormat) {
    ptFormat = tconv.FormatFromExt(testfile);
    if (!ptFormat) {
      cerr << program_name << ": cannot read test file format!" << endl;
      return (-1);
    }
  }

  // Run RMSD!!
  rconv.SetInStream(&rfs);
  tconv.SetInStream(&tfs);
  rconv.SetInFormat(prFormat);
  tconv.SetInFormat(ptFormat);

  cout << "**Starting calcrmsd " << endl;
  cout << "..Reference file = " << reference << endl;
  cout << "..Test file = " << testfile << "\n\n";

  CalcRMSD calcrmsd(rconv, tconv, rmsd_cutoff);
  calcrmsd.Run();

  cout << endl << "**Finishing calcrmsd" << endl;

}

CalcRMSD::CalcRMSD(OBConversion &rconv, OBConversion &tconv, double rmsd_cutoff):
  rconv(rconv), tconv(tconv), rmsd_cutoff(rmsd_cutoff)
{
}

void CalcRMSD::Run()
{

  OBAlign align;
  unsigned int cutoff_passed = 0;
  unsigned int N = 0;
  OBMol rmol, tmol;
  vector<double> rmsd;

  const double binvals_arr[8] = {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 100.0};
  vector<double> binvals(binvals_arr, binvals_arr+8);

  // For every test_conf there is a reference, but not v.v.
  bool more_to_come = tconv.Read(&tmol);
  while (more_to_come) {
    N++;
    string title = tmol.GetTitle();

    // Get the corresponding reference
    rconv.Read(&rmol);
    while (rmol.GetTitle() != title) {
      cout << "..Molecule " << N
           << "\n..title = " << rmol.GetTitle()
           << "\n..number of confs = 0\n";
      N++;
      rconv.Read(&rmol);
    }
    cout << "..Molecule " << N << "\n..title = " << rmol.GetTitle() << "\n";

    align.SetRefMol(rmol);
    rmsd.clear();
    while (tmol.GetTitle() == title) {
      align.SetTargetMol(tmol);
      align.Align();
      rmsd.push_back(align.GetRMSD());
      more_to_come = tconv.Read(&tmol);
      if (!more_to_come) break;
    }

    sort(rmsd.begin(), rmsd.end());

    cout << "..number of confs = " << rmsd.size() << "\n";
    cout << "..minimum rmsd = " << rmsd.at(0) << "\n";

    int bin_idx = 0;
    vector<int> bins(binvals.size());
    for(vector<double>::iterator it=rmsd.begin(); it!=rmsd.end(); ++it) {
      while (*it > binvals[bin_idx])
        bin_idx++;
      bins[bin_idx]++;
    }

    vector<int> cumbins(bins);
    for(int i=1; i<8; ++i)
      cumbins[i] += cumbins[i-1];

    cout << "..confs less than cutoffs: ";
    cout << binvals[0];
    for (int i=1; i < binvals.size(); i++)
      cout << " " << binvals[i];
    cout << "\n";

    cout << ".." << cumbins[0];
    for (int i=1; i < cumbins.size(); i++)
      cout << " " << cumbins[i];
    cout << "\n";

    cout << "..cutoff (" << rmsd_cutoff << ") passed = ";
    if (rmsd.at(0) <= rmsd_cutoff) {
      cout << " Yes\n";
      cutoff_passed++;
    }
    else
      cout << " No\n";
  cout << "\n";

  }

  cout << "\n**Summary\n..number of molecules = " << N
       << "\n..less than cutoff(" << rmsd_cutoff << ") = " << cutoff_passed << "\n";

}

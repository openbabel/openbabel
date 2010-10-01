// Old-style test for Gasteiger charge model

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/chargemodel.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "charge-gasteiger.txt";
  string dipole_file = testdatadir + "dipole-gasteiger.txt";
  string molecules_file = testdatadir + "forcefield.sdf";
#else
  string results_file = "files/charge-gasteiger.txt";
  string dipole_file = "files/dipole-gasteiger.txt";
  string molecules_file = "files/forcefield.sdf";
#endif

void GenerateCharges();

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif

  if (argc != 1)
    {
      if (strncmp(argv[1], "-g", 2))
        {
          cout << "Usage: charge-gasteiger" << endl;
          return 0;
        }
      else
        {
          GenerateCharges();
          return 0;
        }
    }

  cout << "# Testing GASTEIGER Charge Model..." << endl;

  std::ifstream mifs;
  if (!SafeOpen(mifs, molecules_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << molecules_file << endl;
      return -1; // test failed
    }

  std::ifstream rifs;
  if (!SafeOpen(rifs, results_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << results_file << endl;
      return -1; // test failed
    }

  std::ifstream difs;
  if (!SafeOpen(difs, dipole_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << dipole_file << endl;
      return -1; // test failed
    }

  char buffer[BUFF_SIZE];
  vector<string> vs;
  OBMol mol;
  OBConversion conv(&mifs, &cout);
  unsigned int currentTest = 0;
  vector3 dipoleMoment, result;
  
  std::vector<double> partialCharges;

  if(! conv.SetInAndOutFormats("SDF","SDF"))
    {
      cout << "Bail out! SDF format is not loaded" << endl;
      return -1; // test failed
    }
    
  OBChargeModel *pCM = OBChargeModel::FindType("gasteiger");

  if (pCM == NULL) {
    cerr << "Bail out! Cannot load charge model!" << endl;
    return -1; // test failed
  }

  while(mifs)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;
      if (!difs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data" << endl;
          return -1; // test failed
        }
        
      if (!pCM->ComputeCharges(mol)) {
        cout << "Bail out! could not compute charges on " << mol.GetTitle() << endl;
        return -1; // test failed
      }
      partialCharges = pCM->GetPartialCharges();

      // compare the calculated energy to our reference data
      tokenize(vs, buffer);
      if (vs.size() < 3)
        return -1;

      dipoleMoment.SetX(atof(vs[0].c_str()));
      dipoleMoment.SetY(atof(vs[1].c_str()));
      dipoleMoment.SetZ(atof(vs[2].c_str()));
      result = pCM->GetDipoleMoment(mol) - dipoleMoment;
                        
      if ( fabs(result.length_2()) > 1.0e-4)
        {
          cout << "not ok " << ++currentTest << " # calculated dipole incorrect "
               << " for molecule " << mol.GetTitle() << '\n';
        }
      else
        cout << "ok " << ++currentTest << " # dipole\n";

      
      FOR_ATOMS_OF_MOL(atom, mol) {
        if (!rifs.getline(buffer,BUFF_SIZE)) {
          cout << "Bail out! Cannot read reference data\n";
          return -1; // test failed
        }
        
        if ( fabs(atom->GetPartialCharge() - atof(buffer)) > 1.0e-3 ) {
          cout << "not ok " << ++currentTest << " # calculated charge incorrect "
               << " for molecule " << mol.GetTitle() << '\n';
          cout << "# atom " << atom->GetIdx() << " expected " << buffer << " got "
               << atom->GetPartialCharge() << '\n';
        } else {
          cout << "ok " << ++currentTest << " # charge\n";
        }
        
      }

    }

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateCharges()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs, molecules_file.c_str()))
    return;

  std::ofstream rofs;
  if (!SafeOpen(rofs, results_file.c_str()))
    return;

  std::ofstream dofs;
  if (!SafeOpen(dofs, dipole_file.c_str()))
    return;

  OBMol mol;
  OBConversion conv(&ifs, &cout);
  char buffer[BUFF_SIZE];
  
  if(! conv.SetInAndOutFormats("SDF","SDF"))
    {
      cerr << "SDF format is not loaded" << endl;
      return;
    }

  OBChargeModel *pCM = OBChargeModel::FindType("gasteiger");

  if (pCM == NULL) {
    cerr << "Cannot load charge model!" << endl;
    return;
  }

  std::vector<double> partialCharges;
  vector3 dipoleMoment;
  for (;ifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      if (pCM->ComputeCharges(mol)) {
        partialCharges = pCM->GetPartialCharges();
      }

      // write out the dipole moment
      dipoleMoment = pCM->GetDipoleMoment(mol);
      sprintf(buffer, "%15.5f%15.5f%15.5f\n", dipoleMoment.x(), dipoleMoment.y(), dipoleMoment.z());
      dofs << buffer;
      
      // and write all the partial charges
      FOR_ATOMS_OF_MOL(atom, mol) {
        sprintf(buffer, "%15.5f\n", atom->GetPartialCharge());
        rofs << buffer;
      }
    }

	cerr << "Charges written successfully" << endl;
  return;
}

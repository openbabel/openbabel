#include "mol.h"
#include "fingerprint.h"

using namespace std;
using namespace OpenBabel;

int main (int argc, char **argv)
{
  
  if (argc != 3) {
    cerr << "Usage: obtanimoto <query molecule> <fingerprints>\n";
    exit (1);
  } 
  char *program_name = argv[0];
  char *FileIn = argv[1], *FptIn = argv[2];
  io_type inFileType = UNDEFINED, outFileType = FINGERPRINT;
  ifstream ifs;
  unsigned int i;
  // Find Input filetype
  if (extab.CanReadExtension(FileIn))
    inFileType = extab.FilenameToType(FileIn);
  else
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }
  //if (extab.CanWriteExtension(FileIn))
  //  outFileType = extab.FilenameToType(FileIn);
  //else
  //  {
  //    cerr << program_name << ": cannot write input format!" << endl;
  //    exit (-1);
  //  }

  
  OBMol mol(inFileType, outFileType);
  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input molecule file!" << endl;
      exit (-1);
    }

  ifs >> mol;

  ifs.close();
  // hash the molecule fragments
  fingerprint fpt(mol.GetTitle());
  fpt.HashMol(mol);

  ifs.open(FptIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input fingerprint file!" << endl;
      exit (-1);
    }
  string line, name;
  OBBitVec fingerprint(1024);
  double t;
  while ( 1 ) {
    if ( ! getline(ifs,line,'\n') )
      break;
    name= line;
    //erase '>'
    name.erase(0,1);
    // Get the fingerprint line
    if ( !getline(ifs,line,'\n')) {
      cerr << "Error: Unable to read second fingerprint line\n.";
      return 1;
    }
    if (line.size() != 1024) {
      cerr << "Error: bad fingerprint size: " << line.size();
      return 1;
    }
    // Set the fingerprint bits
    fingerprint.Clear();
    for(i=0; i<line.size(); i++) {
      if (line[i] == '1')
	fingerprint.SetBitOn(i);
    }
    t = Tanimoto(fpt.GetFpt(),fingerprint);
    cout << name << " " << t << endl;
  }
  return 0;
}

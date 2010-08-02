#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>
#include <iomanip>

using namespace OpenBabel;
int main(int argc,char **argv)
{
  if(argc==2)
  {
    //Make a fingerprint of type FP2(the default)
    //Equivalent to babel filename.xxx -ofpt -xo

    const char* filename = argv[1];
    OBConversion conv;

    //Read in a molecule from a file of any format
    OBMol mol;
    OBFormat* pFormat = conv.FormatFromExt(filename);
    if(pFormat && conv.SetInFormat(pFormat) && conv.ReadFile(&mol, filename))
    {
      //Get a pointer to the OBFingerprint class with id=="FP2"
      OBFingerprint* fptype = OBFingerprint::FindType("FP2");

      //Make a fingerprint of the default size from the molecule
      std::vector<unsigned> fptvec;
      if(fptype->GetFingerprint(&mol, fptvec))
      {
        //Output the fingerprint as hexadecimal, with six 32bit words per line
        for(int i=fptvec.size()-1;i>=0;i--)
        {
          std::cout << std::hex << std::setfill('0') << std::setw(8) << fptvec[i] << " " ;
          if((fptvec.size()-i)%6==0)
            std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }
  }
  else
    std::cout << "Usage:\n" << argv[0] << " filename.xxx" << std::endl;

  return 0;
}

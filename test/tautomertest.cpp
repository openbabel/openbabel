#include "obtest.h"
#include <openbabel/tautomer.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <algorithm>

using namespace OpenBabel;

/**
 * Check the number of enumerated tautomers.
 */
void testEnumerateTautomers(const std::string &smiles, int numTautomers)
{
  class Functor : public UniqueTautomerFunctor
  {
    public:
      int numTautomers;

      Functor() : numTautomers(0) {}
      void operator()(OBMol*, const std::string&)
      {
        numTautomers++;
      }
  };

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smiles);

  Functor functor;
  EnumerateTautomers(&mol, functor);

  OB_COMPARE( functor.numTautomers, numTautomers );
}

/**
 * Check that each enumerated tautomer has the same canonical tautomer.
 */
void testCanonicalTautomers(const std::string &smiles)
{
  class Functor : public TautomerFunctor
  {
    public:
      std::vector<std::string> tautomers;

      void operator()(OBMol *mol)
      {
        OBConversion conv;
        conv.SetOutFormat("can");
        tautomers.push_back(conv.WriteString(mol, true));
      }
  };

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");
  conv.ReadString(&mol, smiles);

  // enumerate all tautomers
  Functor functor;
  EnumerateTautomers(&mol, functor);

  // check to make sure all tautomers result in the same canonical tautomer
  const std::vector<std::string> &tautomers = functor.tautomers;
  if (tautomers.empty())
    return;

  std::vector<std::string> canonicalTautomers;
  for (std::size_t i = 0; i < tautomers.size(); ++i) {
    OBMol mol2;
    conv.ReadString(&mol2, tautomers[i]);
    CanonicalTautomer(&mol2);
    canonicalTautomers.push_back(conv.WriteString(&mol2, true));
  }

  canonicalTautomers.erase(std::unique(canonicalTautomers.begin(), canonicalTautomers.end()), canonicalTautomers.end());

  OB_COMPARE(canonicalTautomers.size(), unsigned(1));
}

/**
 * Verify that the canonical tautomer is the same as the one provided.
 */
void testVerifyCanonicalTautomer(const std::string &smiles, const std::string &expected)
{
  // Read the SMILES string
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");
  conv.ReadString(&mol, smiles);

  // Get the canonical tautomer
  CanonicalTautomer(&mol);

  std::string experimental = conv.WriteString(&mol, true);

  OB_COMPARE(expected, experimental);
}


int tautomertest(int argc, char* argv[]) {
  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif


    int defaultchoice = 1;

    int choice = defaultchoice;

    if (argc > 1) {
      if(sscanf(argv[1], "%d", &choice) != 1) {
        std::cout << "Couldn't parse that input as a number\n";
        return -1;
      }
    }


    // guanine
    switch (choice) {
    case 1:
      // guanine
      testEnumerateTautomers("Nc1nc2ncnc2c([nH]1)O", 15);
      break;
    case 2:
      // each output for guanine, should also have 15 tautomers
      testEnumerateTautomers("Nc1nc(=O)c2=c([nH]1)[nH]cn2", 15);
      testEnumerateTautomers("Nc1nc(=O)c2=c([nH]1)nc[nH]2", 15);
      testEnumerateTautomers("Oc1nc(N)[nH]c2-c1ncn2", 15);
      testEnumerateTautomers("Nc1[nH]c(=O)c2=c(n1)[nH]cn2", 15);
      testEnumerateTautomers("Nc1nc(O)c2-c(n1)[nH]cn2", 15);
      testEnumerateTautomers("Nc1nc2=c(c(=O)[nH]1)[nH]cn2", 15);
      testEnumerateTautomers("Nc1nc(O)c2-c(n1)nc[nH]2", 15);
      testEnumerateTautomers("Nc1nc2-c(c([nH]1)O)ncn2", 15);
      testEnumerateTautomers("N=c1[nH]c(=O)c2=c([nH]1)[nH]cn2", 15);
      testEnumerateTautomers("N=c1nc(O)c2=c([nH]1)[nH]cn2", 15);
      testEnumerateTautomers("N=c1[nH]c(=O)c2=c([nH]1)nc[nH]2", 15);
      testEnumerateTautomers("N=c1nc(O)c2=c([nH]1)nc[nH]2", 15);
      testEnumerateTautomers("N=c1[nH]c(O)c2-c([nH]1)ncn2", 15);
      testEnumerateTautomers("N=c1[nH]c(O)c2-c(n1)[nH]cn2", 15);
      testEnumerateTautomers("N=c1nc2-c(c([nH]1)O)[nH]cn2", 15);
      break;
    case 3:
      // each output for guanine should have the same canonical tautomer
      testCanonicalTautomers("Nc1nc2ncnc2c([nH]1)O");
      testCanonicalTautomers("Nc1nc(=O)c2=c([nH]1)[nH]cn2");
      testCanonicalTautomers("Nc1nc(=O)c2=c([nH]1)nc[nH]2");
      testCanonicalTautomers("Oc1nc(N)[nH]c2-c1ncn2");
      testCanonicalTautomers("Nc1[nH]c(=O)c2=c(n1)[nH]cn2");
      testCanonicalTautomers("Nc1nc(O)c2-c(n1)[nH]cn2");
      testCanonicalTautomers("Nc1nc2=c(c(=O)[nH]1)[nH]cn2");
      testCanonicalTautomers("Nc1nc(O)c2-c(n1)nc[nH]2");
      testCanonicalTautomers("Nc1nc2-c(c([nH]1)O)ncn2");
      testCanonicalTautomers("N=c1[nH]c(=O)c2=c([nH]1)[nH]cn2");
      testCanonicalTautomers("N=c1nc(O)c2=c([nH]1)[nH]cn2");
      testCanonicalTautomers("N=c1[nH]c(=O)c2=c([nH]1)nc[nH]2");
      testCanonicalTautomers("N=c1nc(O)c2=c([nH]1)nc[nH]2");
      testCanonicalTautomers("N=c1[nH]c(O)c2-c([nH]1)ncn2");
      testCanonicalTautomers("N=c1[nH]c(O)c2-c(n1)[nH]cn2");
      testCanonicalTautomers("N=c1nc2-c(c([nH]1)O)[nH]cn2");
      break;
    case 4:
      // benzene
      testEnumerateTautomers("c1ccccc1", 1);
      testEnumerateTautomers("C1=CC=CC=C1", 1);
      testEnumerateTautomers("C=1C=CC=CC1", 1);
      testVerifyCanonicalTautomer("c1ccccc1", "c1ccccc1");
      testVerifyCanonicalTautomer("C1=CC=CC=C1", "c1ccccc1");
      testVerifyCanonicalTautomer("C=1C=CC=CC1", "c1ccccc1");
      break;
    case 5:
      // vinyl ether
      testEnumerateTautomers("C=COC=C", 1);
      testVerifyCanonicalTautomer("C=COC=C", "C=COC=C");
      break;
    case 6:
      // 6-mercaptopurine
      testEnumerateTautomers("C1=NC2=C(N1)C(=S)N=CN2", 8);
      testCanonicalTautomers("C1=NC2=C(N1)C(=S)N=CN2");
      testVerifyCanonicalTautomer("C1=NC2=C(N1)C(=S)N=CN2", "Sc1[nH]cnc2-c1ncn2");
      break;
    case 7:
      // allopurinol
      testEnumerateTautomers("C1=C2C(=NC=NC2=O)NN1", 9);
      testCanonicalTautomers("C1=C2C(=NC=NC2=O)NN1");
      testVerifyCanonicalTautomer("C1=C2C(=NC=NC2=O)NN1", "Oc1[nH]cnc2-c1cnn2");
      break;
    case 8:
      // chlorzoxazone
      testEnumerateTautomers("C1=CC2=C(C=C1Cl)NC(=O)O2", 2);
      testCanonicalTautomers("C1=CC2=C(C=C1Cl)NC(=O)O2");
      testVerifyCanonicalTautomer("C1=CC2=C(C=C1Cl)NC(=O)O2", "Clc1ccc2c(c1)[nH]c(=O)o2");
      break;
    case 9:
      // 2-mercaptobenzothiazole
      testEnumerateTautomers("C1=CC=C2C(=C1)NC(=S)S2", 2);
      testCanonicalTautomers("C1=CC=C2C(=C1)NC(=S)S2");
      testVerifyCanonicalTautomer("C1=CC=C2C(=C1)NC(=S)S2", "Sc1nc2c(s1)cccc2");
      break;
    case 10:
      // pemoline
      testEnumerateTautomers("C1=CC=C(C=C1)C2C(=O)N=C(O2)N", 3);
      testCanonicalTautomers("C1=CC=C(C=C1)C2C(=O)N=C(O2)N");
      testVerifyCanonicalTautomer("C1=CC=C(C=C1)C2C(=O)N=C(O2)N", "O=C1N=C(OC1c1ccccc1)N");
      break;
    case 11:
      // amitrole
      testEnumerateTautomers("C1=NNC(=N1)N", 5);
      testCanonicalTautomers("C1=NNC(=N1)N");
      testVerifyCanonicalTautomer("C1=NNC(=N1)N", "Nc1ncn[nH]1");
      break;
    case 12:
      // purpald
      testEnumerateTautomers("C1(=S)NN=C(N1N)NN", 4);
      testCanonicalTautomers("C1(=S)NN=C(N1N)NN");
      testVerifyCanonicalTautomer("C1(=S)NN=C(N1N)NN", "Nn1c(NN)nnc1S");
      break;
    case 13:
      // thiotetronic acid
      testEnumerateTautomers("O=C1CSC(=O)C1", 1);
      testCanonicalTautomers("O=C1CSC(=O)C1");
      testVerifyCanonicalTautomer("O=C1CSC(=O)C1", "O=C1SCC(=O)C1");
      break;
    case 14:
      // pyrithione
      testEnumerateTautomers("C1=CC(=S)N(C=C1)O", 1);
      testCanonicalTautomers("C1=CC(=S)N(C=C1)O");
      testVerifyCanonicalTautomer("C1=CC(=S)N(C=C1)O", "S=c1ccccn1O");
      break;
    case 15:
      // iodothiouracil
      testEnumerateTautomers("C1=C(C(=O)NC(=S)N1)I", 6);
      testCanonicalTautomers("C1=C(C(=O)NC(=S)N1)I");
      testVerifyCanonicalTautomer("C1=C(C(=O)NC(=S)N1)I", "O=c1nc(S)[nH]cc1I");
      break;
    case 16:
      // divicine
      testEnumerateTautomers("n1c(N)nc(O)c(O)c1(N)", 9);
      testCanonicalTautomers("n1c(N)nc(O)c(O)c1(N)");
      testVerifyCanonicalTautomer("n1c(N)nc(O)c(O)c1(N)", "O=c1nc(N)[nH]c(c1O)N");
      break;
    case 17:
      // 2-thiouracil
      testEnumerateTautomers("C1=CNC(=S)NC1=O", 6);
      testCanonicalTautomers("C1=CNC(=S)NC1=O");
      testVerifyCanonicalTautomer("C1=CNC(=S)NC1=O", "Oc1ccnc(=S)[nH]1");
      break;
    case 18:
      // flucytosine
      testEnumerateTautomers("C1=NC(=O)NC(=C1F)N", 6);
      testCanonicalTautomers("C1=NC(=O)NC(=C1F)N");
      testVerifyCanonicalTautomer("C1=NC(=O)NC(=C1F)N", "N=c1nc(O)[nH]cc1F");
      break;
    case 19:
      // citrazinic acid
      testEnumerateTautomers("C1=C(C=C(NC1=O)O)C(=O)O", 2);
      testCanonicalTautomers("C1=C(C=C(NC1=O)O)C(=O)O");
      testVerifyCanonicalTautomer("C1=C(C=C(NC1=O)O)C(=O)O", "OC(=O)c1cc(O)[nH]c(=O)c1");
      break;
    case 20:
      // guanoxabenz
      testEnumerateTautomers("C1=CC(=C(C(=C1)Cl)C=NN=C(N)NO)Cl", 3);
      testCanonicalTautomers("C1=CC(=C(C(=C1)Cl)C=NN=C(N)NO)Cl");
      testVerifyCanonicalTautomer("C1=CC(=C(C(=C1)Cl)C=NN=C(N)NO)Cl", "ONC(=N)NN=Cc1c(Cl)cccc1Cl");
      break;
    case 21:
      // tenoxicam
      testEnumerateTautomers("CN1/C(=C(\\NC2=CC=CC=N2)/O)/C(=O)C3=C(S1(=O)=O)C=CS3", 5);
      testCanonicalTautomers("CN1/C(=C(\\NC2=CC=CC=N2)/O)/C(=O)C3=C(S1(=O)=O)C=CS3");
      testVerifyCanonicalTautomer("CN1/C(=C(\\NC2=CC=CC=N2)/O)/C(=O)C3=C(S1(=O)=O)C=CS3", "OC(=C1C(=O)c2sccc2S(=O)(=O)N1C)Nc1ccccn1");
      break;
    case 22:
      // mitoguazone
      testEnumerateTautomers("CC(=NN=C(N)N)C=NN=C(N)N", 8);
      testCanonicalTautomers("CC(=NN=C(N)N)C=NN=C(N)N");
      testVerifyCanonicalTautomer("CC(=NN=C(N)N)C=NN=C(N)N", "CC(=CN=NC(=N)N)NNC(=N)N");
      break;
    case 23:
      // leucopterin
      testEnumerateTautomers("C12=C(NC(=O)C(=O)N1)N=C(NC2=O)N", 33);
      testCanonicalTautomers("C12=C(NC(=O)C(=O)N1)N=C(NC2=O)N");
      testVerifyCanonicalTautomer("C12=C(NC(=O)C(=O)N1)N=C(NC2=O)N", "Nc1nc(=O)c2c([nH]1)[nH]c(=O)c(=O)[nH]2");
      break;
    case 24:
      // methimazole
      testEnumerateTautomers("CN1C=CNC1=S", 2);
      testCanonicalTautomers("CN1C=CNC1=S");
      testVerifyCanonicalTautomer("CN1C=CNC1=S", "Cn1ccnc1S");
      break;
    case 25:
      // ciclopirox
      testEnumerateTautomers("CC1=CC(=O)N(C(=C1)C2CCCCC2)O", 1);
      testCanonicalTautomers("CC1=CC(=O)N(C(=C1)C2CCCCC2)O");
      testVerifyCanonicalTautomer("CC1=CC(=O)N(C(=C1)C2CCCCC2)O", "Cc1cc(C2CCCCC2)n(c(=O)c1)O");
      break;
    case 26:
      // violuric acid
      testEnumerateTautomers("C1(=C(NC(=O)NC1=O)O)N=O", 10);
      testCanonicalTautomers("C1(=C(NC(=O)NC1=O)O)N=O");
      testVerifyCanonicalTautomer("C1(=C(NC(=O)NC1=O)O)N=O", "ON=C1C(=NC(=O)N=C1O)O");
      break;
    case 27:
      // methisazone
      testEnumerateTautomers("CN1C2=CC=CC=C2C(=NNC(=S)N)C1=O", 5);
      testCanonicalTautomers("CN1C2=CC=CC=C2C(=NNC(=S)N)C1=O");
      testVerifyCanonicalTautomer("CN1C2=CC=CC=C2C(=NNC(=S)N)C1=O", "NC(=S)NN=C1c2ccccc2N(C1=O)C");
      break;
    case 28:
      // antralin
      testEnumerateTautomers("Oc1cccc2c1C(=O)c1c(C2)cccc1O", 2);
      testCanonicalTautomers("Oc1cccc2c1C(=O)c1c(C2)cccc1O");
      testVerifyCanonicalTautomer("Oc1cccc2c1C(=O)c1c(C2)cccc1O", "O=C1C=CC=C2C1=C(O)c1c(C2)cccc1O");
      break;
    default:
      std::cout << "Test number " << choice << " does not exist!\n";
      return -1;
    }

  return 0;
}

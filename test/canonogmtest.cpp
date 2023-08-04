#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>


#include <iostream>
#include <vector>
#include <algorithm>
#include "../src/formats/smilesformat.cpp"

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
    string ogmtestdatadir = TESTDATADIR;
    //string complete_dataset_file = ogmtestdatadir + "ogm/complete_dataset.smi"; //"C:/TFG/dataset/complete_dataset.smi";
    //string dataset_file =         ogmtestdatadir + "ogm/dataset_to_test_small.smi"; //"C:/TFG/dataset/output/test/files/dataset_to_test_small.smi";

    //string dataset_file = ogmtestdatadir + "ogm/prueba.smi";

    string dataset_file =           ogmtestdatadir + "ogm/dataset_to_test.smi";
    string canon_file =             ogmtestdatadir + "ogm/canon_output.smi";
    string random_file =            ogmtestdatadir + "ogm/random_to_test.smi";
    string random_canon_file =      ogmtestdatadir + "ogm/random_canon_output.smi";
    string canon_over_canon_file =  ogmtestdatadir + "ogm/canon_over_canon_output.smi";
    string metal_selection_file =   ogmtestdatadir + "ogm/metal_selection.smi";
#endif



//Test unitario para comprobar la seleccion del primer metal
int SelectCanonMetal() {
    cout << endl << "# Testing Metal Selection in Canonical Algorithm...  \n";

    //Abrimos el fichero con los smiles a testear
    std::ifstream mifs;
    if (!SafeOpen(mifs, metal_selection_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << metal_selection_file << endl;
        return -1; // test failed
    }

    OutOptions options(true, false,
        false,
        false, false,
        nullptr);

    OBMol2Cansmi m2s(options);

    std::string buffer, smiles;
    std::vector<int> idx;
    std::string delimiter = "?";
    while (mifs.good())
    {
        if (std::getline(mifs, buffer)) { //Sacamos la linea entera y luego la parseamos

            if ((buffer[0] == '#' && buffer[1] == '#') || buffer.empty()) //Es una linea de comentarios o linea vacia (probablemente al final del fichero)
                continue;
            else {
                //Parse string to get smiles and metal_idx
                size_t pos = 0;
                std::string token;
                if (pos = buffer.find(delimiter)) {
                    smiles = buffer.substr(0, pos);
                    buffer.erase(0, pos + delimiter.length());
                }

                while ((pos = buffer.find(delimiter)) != std::string::npos) {
                    token = buffer.substr(0, pos);
                    idx.push_back(std::stoi(token));
                    buffer.erase(0, pos + delimiter.length());
                }
            }

            //Safety checks for invalids idx
            bool toContinue = false;
            for (auto c : idx) {
                if (c <= 0)
                    toContinue = true;
            }
            if (toContinue) {
                idx.clear();
                buffer.clear();
                smiles.clear();
                continue;
            }

            //Ahora ya hacemos el test de seleccion para esta molecula parseada
            //...
            bool goodTest = false;
            OBConversion conv;
            conv.SetInFormat("smi");
            conv.SetOutFormat("smi");
            OBMol mol;
            conv.ReadString(&mol, smiles);


            m2s.Init(&mol, false, &conv);

            OBAtom* _startatom = m2s.SelectRootAtomOgm(mol, &conv);

            cout << "SMILES: " << smiles;
            cout << "\nAtomo(s) correcto(s):";
            for (auto c : idx) {
                std::cout << " " << c;
            }
            cout << " \nAtomo escogido : " << _startatom->GetIdx() << "\n\n";

            goodTest = (std::find(idx.begin(), idx.end(), _startatom->GetIdx()) != idx.end());
            OB_REQUIRE(goodTest);


            //Limpiamos las variables
            idx.clear();
            buffer.clear();
            smiles.clear();
            mol.Clear();
            goodTest = false;
        }

    }

    mifs.close();
    mifs.clear();



    return 0;
}


//Test funcional para comprobar la robustez del canonizado
int RandomCanonStandardLabels() {
    cout << endl << "# Testing Canon Persistance when using random SMILES with Standard Labels algorithm...\n";


    //Convertimos los SMILES originales a random SMILES usando la opcion -C anticanonical
    OBConversion convAntiC;
    convAntiC.SetInAndOutFormats("smi", "smi");
    convAntiC.AddOption("C"); //Anticanonical smiles
    vector<string> FileList, OutputFileList;
    string OutputFileName;
    FileList.push_back(dataset_file);
    OutputFileName = random_file;
    int count = convAntiC.FullConvert(FileList, OutputFileName, OutputFileList);
    OB_REQUIRE((count > 0));


    //Convertimos los random SMILES en canonicos
    OBConversion convCanon;
    convCanon.SetInAndOutFormats("smi", "smi");
    vector<string> FileList2, OutputFileList2;
    string OutputFileName2;
    FileList2.push_back(random_file);
    OutputFileName2 = random_canon_file;
    count = convCanon.FullConvert(FileList2, OutputFileName2, OutputFileList2);
    OB_REQUIRE((count > 0));

    //Convertimos los original SMILES en canonicos
    OBConversion convCanon2;
    convCanon2.SetInAndOutFormats("smi", "smi");
    vector<string> FileList3, OutputFileList3;
    string OutputFileName3;
    FileList3.push_back(dataset_file);
    OutputFileName3 = canon_file;
    count = convCanon2.FullConvert(FileList3, OutputFileName3, OutputFileList3);
    OB_REQUIRE((count > 0));

    std::ifstream mifs, mifs2, mifs3, mifs4;
    if (!SafeOpen(mifs, dataset_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << dataset_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs2, random_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << random_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs3, canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << canon_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs4, random_canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << random_canon_file << endl;
        return -1; // test failed
    }


    //Leemos ambos ficheros y comparamos los resultados
    unsigned int currentTest = 0;
    std::string canon, random_canon, random, original;
    while (mifs.good() && mifs2.good() && mifs3.good())
    {
        if (std::getline(mifs, original) && std::getline(mifs2, random) && std::getline(mifs3, canon) && std::getline(mifs4, random_canon)) {
            cout << "Original SMILES for \t\t\t" << ++currentTest << ": " << original << "\n";
            cout << "Random SMILES for \t\t\t" << currentTest << ": " << random << "\n";
            cout << "Canon for Original SMILES for \t" << currentTest << ": " << canon << "\n";
            cout << "Canon for Random SMILES for \t" << currentTest << ": " << random_canon << "\n\n";
            OB_COMPARE(random_canon, canon);
        }

    }

    mifs.close();
    mifs.clear();

    mifs2.close();
    mifs2.clear();

    mifs3.close();
    mifs3.clear();

    mifs4.close();
    mifs4.clear();

    return 0;
}

//Test funcional para comprobar la robustez del canonizado
int RandomCanonCanonicalLabels() {
    cout << endl << "# Testing Canon Persistance when using random SMILES with Canonical Labels algorithm...\n";


    //Convertimos los SMILES originales a random SMILES usando la opcion -C anticanonical
    OBConversion convAntiC;
    convAntiC.SetInAndOutFormats("smi", "smi");
    convAntiC.AddOption("C"); //Anticanonical smiles
    vector<string> FileList, OutputFileList;
    string OutputFileName;
    FileList.push_back(dataset_file);
    OutputFileName = random_file;
    int count = convAntiC.FullConvert(FileList, OutputFileName, OutputFileList);
    OB_REQUIRE((count > 0));


    //Convertimos los random SMILES en canonicos
    OBConversion convCanon;
    convCanon.SetInAndOutFormats("smi", "smi");
    convCanon.AddOption("c");
    vector<string> FileList2, OutputFileList2;
    string OutputFileName2;
    FileList2.push_back(random_file);
    OutputFileName2 = random_canon_file;
    count = convCanon.FullConvert(FileList2, OutputFileName2, OutputFileList2);
    OB_REQUIRE((count > 0));

    //Convertimos los original SMILES en canonicos
    OBConversion convCanon2;
    convCanon2.SetInAndOutFormats("smi", "smi");
    convCanon2.AddOption("c");
    vector<string> FileList3, OutputFileList3;
    string OutputFileName3;
    FileList3.push_back(dataset_file);
    OutputFileName3 = canon_file;
    count = convCanon2.FullConvert(FileList3, OutputFileName3, OutputFileList3);
    OB_REQUIRE((count > 0));

    std::ifstream mifs, mifs2, mifs3, mifs4;
    if (!SafeOpen(mifs, dataset_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << dataset_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs2, random_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << random_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs3, canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << canon_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs4, random_canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << random_canon_file << endl;
        return -1; // test failed
    }


    //Leemos ambos ficheros y comparamos los resultados
    unsigned int currentTest = 0;
    std::string canon, random_canon, random, original;
    while (mifs.good() && mifs2.good() && mifs3.good())
    {
        if (std::getline(mifs, original) && std::getline(mifs2, random) && std::getline(mifs3, canon) && std::getline(mifs4, random_canon)) {
            cout << "Original SMILES for \t\t" << ++currentTest << ": " << original << "\n";
            cout << "Random SMILES for \t\t" << currentTest << ": " << random << "\n";
            cout << "Canon for Original SMILES for \t" << currentTest << ": " << canon << "\n";
            cout << "Canon for Random SMILES for \t" << currentTest << ": " << random_canon << "\n\n";
            OB_COMPARE(random_canon, canon);
        }

    }

    mifs.close();
    mifs.clear();

    mifs2.close();
    mifs2.clear();

    mifs3.close();
    mifs3.clear();

    mifs4.close();
    mifs4.clear();

    return 0;
}

//Test funcional para comprobar la robustez del canonizado
int CanonOverCanon() {
    cout << endl << "# Testing Canon Persistance when using a previously canonized SMILES (i.e. aplying 2 times the algorithm)...\n";


    //Convertimos los SMILES de entrada en canonicos
    int count = 0;
    OBConversion convCanon;
    convCanon.SetInAndOutFormats("smi", "smi");
    convCanon.AddOption("c");
    vector<string> FileList, OutputFileList;
    string OutputFileName;
    FileList.push_back(dataset_file);
    OutputFileName = canon_file;
    count = convCanon.FullConvert(FileList, OutputFileName, OutputFileList);
    OB_REQUIRE((count > 0));

    //Convertimos los SMILES resultado de la conversion anterior de nuevo en canonicos
    OBConversion convCanon2;
    convCanon2.SetInAndOutFormats("smi", "smi");
    convCanon2.AddOption("c");
    vector<string> FileList2, OutputFileList2;
    string OutputFileName2;
    FileList2.push_back(canon_file);
    OutputFileName2 = canon_over_canon_file;
    count = convCanon2.FullConvert(FileList2, OutputFileName2, OutputFileList2);
    OB_REQUIRE((count > 0));

    std::ifstream mifs, mifs2, mifs3;
    if (!SafeOpen(mifs, dataset_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << dataset_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs2, canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << canon_file << endl;
        return -1; // test failed
    }
    if (!SafeOpen(mifs3, canon_over_canon_file.c_str()))
    {
        cout << "Bail out! Cannot read file " << canon_over_canon_file << endl;
        return -1; // test failed
    }


    //Leemos ambos ficheros y comparamos los resultados
    unsigned int currentTest = 0;
    std::string canon, canon_over_canon, original;
    while (mifs.good() && mifs2.good() && mifs3.good())
    {
        if (std::getline(mifs, original) && std::getline(mifs2, canon) && std::getline(mifs3, canon_over_canon)) {
            cout << "Original SMILES for \t\t\t" << ++currentTest << ": " << original << "\n";
            cout << "Canon for Original SMILES for \t" << currentTest << ": " << canon << "\n";
            cout << "Canon for Canon SMILES for \t\t" << currentTest << ": " << canon_over_canon << "\n\n";
            OB_COMPARE(canon, canon_over_canon);
        }

    }

    mifs.close();
    mifs.clear();

    mifs2.close();
    mifs2.clear();

    mifs3.close();
    mifs3.clear();

    return 0;
}

int DrawDoubleCpTest() {
    cout << endl << "# Testing Detection and Drawing Double Cp in SVG Depiction...  \n";

    OBConversion conv;
    conv.SetInFormat("smi");
    conv.SetOutFormat("svg");
    OBMol mol;

    vector<string> FileList, OutputFileList;
    string OutputFileName;
    FileList.push_back("-:[Cl-][Au+][P](C=1C=CC=CC1)(C=2C=CC=CC2)[C-]34[CH]5=[CH]6[CH]7=[CH]3[Fe+2]6789%10%1154[CH]=%12[CH]%11=[CH]%10[C-]9([CH]%128)[P]([Au+][Cl-])(C=%13C=CC=CC%13)C=%14C=CC=CC%14");
    OutputFileName = ogmtestdatadir + "ogm/ogm_test_5_dobleCp.svg";

    int count = conv.FullConvert(FileList, OutputFileName, OutputFileList);

    OB_REQUIRE((count > 0));
    return 0;
}






int canonogmtest(int argc, char* argv[]) {

    int defaultchoice = 1;

    int choice = defaultchoice;

    if (argc > 1) {
        if (sscanf(argv[1], "%d", &choice) != 1) {
            printf("Couldn't parse that input as a number\n");
            return -1;
        }
    }
    // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif  
    int result = 0;
    switch (choice) {
    case 1:
        result = SelectCanonMetal();
        break;
    case 2:
        result = RandomCanonStandardLabels();
        break;
    case 3:
        result = RandomCanonCanonicalLabels();
        break;
    case 4:
        result = CanonOverCanon();
        break;
    case 5:
        result = DrawDoubleCpTest();
        break;

        //case N:
        //  YOUR_TEST_HERE();
        //  Remember to update CMakeLists.txt with the number of your test
        //  break;
    default:
        cout << "Test number " << choice << " does not exist!\n";
        return -1;
    }

    return result;
}


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>


extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    using namespace std;
    using namespace OpenBabel;
    OBConversion obconv;
    OpenBabel::OBMol obmol;
    std::string str (reinterpret_cast<const char*>(Data), Size);

    //FUZZ_INPUT_FORMAT is defined at compile time
    if(!obconv.SetInFormat(FUZZ_INPUT_FORMAT)){
        abort();
    }
    obconv.ReadString(&obmol, str);
    return 0;
}


#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <cstdlib>
#include <stdio.h>
#include <sstream>

#include <fuzzer/FuzzedDataProvider.h>

std::vector<std::string> getAllFormats() {
    std::vector<std::string> formats;
    OpenBabel::OBPlugin::ListAsVector("formats", "ids", formats);
    std::sort(formats.begin(), formats.end());
    return formats;
}

const char *randomFormat(FuzzedDataProvider &fdp) {
    static std::vector<std::string> formats = getAllFormats();
    return formats[fdp.ConsumeIntegralInRange<int>(0, formats.size() - 1)].c_str();
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    using namespace OpenBabel;
    obErrorLog.StopLogging();

    FuzzedDataProvider fdp(Data, Size);
    OBConversion conv;
    conv.SetInFormat(randomFormat(fdp));
    conv.SetOutFormat(randomFormat(fdp));

    std::ostringstream outStream;
    std::istringstream inStream(fdp.ConsumeRandomLengthString());

    try {
        conv.Convert(&inStream, &outStream);
    } catch (...) {
        // no error handling
    }

    return 0;
}

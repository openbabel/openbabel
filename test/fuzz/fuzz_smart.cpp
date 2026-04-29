#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <fuzzer/FuzzedDataProvider.h>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    if (Size < 2) return 0;
    
    FuzzedDataProvider fdp(Data, Size);
    
    // Use part of the data for SMARTS pattern
    std::string smarts_pattern = fdp.ConsumeRandomLengthString();
    
    // Use the rest for a SMILES molecule
    std::string smiles = fdp.ConsumeRemainingBytesAsString();
    
    if (smarts_pattern.empty() || smiles.empty()) return 0;

    OpenBabel::obErrorLog.StopLogging();

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;
    if (conv.SetInFormat("SMILES")) {
        if (conv.ReadString(&mol, smiles)) {
            OpenBabel::OBSmartsPattern sp;
            if (sp.Init(smarts_pattern)) {
                sp.Match(mol);
            }
        }
    }

    return 0;
}

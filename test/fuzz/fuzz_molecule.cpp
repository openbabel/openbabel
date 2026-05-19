#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

#include <vector>
#include <string>

// Vendored alongside this file (LLVM Apache-2.0 WITH LLVM-exception)
// so the harness compiles outside of an LLVM/libfuzzer toolchain — in
// particular under GCC, where the libfuzzer headers aren't shipped.
#include "FuzzedDataProvider.h"

using namespace OpenBabel;

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    if (Size < 2) return 0;
    
    FuzzedDataProvider fdp(Data, Size);
    
    // List of common input formats in OpenBabel
    static const std::vector<std::string> formats = {
        "smi", "smiles", "sdf", "mol", "pdb", "ent", "xyz", "cml", "inchi", "mna", "can"
    };
    
    int format_idx = fdp.ConsumeIntegralInRange<int>(0, formats.size() - 1);
    std::string format = formats[format_idx];
    std::string input = fdp.ConsumeRandomLengthString();
    
    if (input.empty()) return 0;

    obErrorLog.StopLogging();

    OBMol mol;
    OBConversion conv;
    if (conv.SetInFormat(format.c_str())) {
        if (conv.ReadString(&mol, input)) {
            // Randomly perform various operations
            
            if (fdp.ConsumeBool()) {
                // Perception routines
                mol.FindSSSR();
                mol.FindLSSR();
                mol.FindRingAtomsAndBonds();
                mol.PerceiveBondOrders();
                mol.FindAngles();
                mol.FindTorsions();
                mol.GetFormula();
                mol.GetMolWt();
                mol.GetExactMass();
            }
            
            if (fdp.ConsumeBool()) {
                // Hydrogens
                if (fdp.ConsumeBool()) {
                    mol.AddHydrogens();
                } else {
                    mol.DeleteHydrogens();
                }
            }

            if (fdp.ConsumeBool()) {
                // 3D Builder
                OBBuilder builder;
                builder.Build(mol);
            }
            
            if (fdp.ConsumeBool()) {
                // Force Fields
                static const std::vector<std::string> ff_names = {
                    "MMFF94", "UFF", "GAFF", "Ghemical"
                };
                int ff_idx = fdp.ConsumeIntegralInRange<int>(0, ff_names.size() - 1);
                std::string ff_name = ff_names[ff_idx];
                OBForceField *ff = OBForceField::FindType(ff_name.c_str());
                if (ff && ff->Setup(mol)) {
                    ff->Energy();
                    uint8_t min_steps = fdp.ConsumeIntegralInRange<uint8_t>(0, 10);
                    if (min_steps > 0) {
                        ff->SteepestDescent(min_steps);
                        ff->ConjugateGradients(min_steps);
                    }
                }
            }
        }
    }

    return 0;
}

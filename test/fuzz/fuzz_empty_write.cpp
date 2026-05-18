// Exercise every writable format with an empty (zero-atom) molecule.
//
// Zero-atom molecules are a recurring source of bugs in writers (and
// some readers that allocate-then-bail): missing null checks, leaks on
// early-return paths, off-by-ones when iterating empty arrays, etc.
// This harness enumerates every writable OBFormat once per invocation
// and asks it to serialise an empty OBMol, so one corpus seed touches
// every writer. The fuzz bytes additionally drive optional title and
// OBPairData content, letting libfuzzer probe writers that emit those.
//
// When LIB_FUZZING_ENGINE is unset, standalone_main.cc supplies main()
// and this also runs as a plain CTest regression case.
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>

#include <algorithm>
#include <string>
#include <vector>

// Vendored alongside this file (LLVM Apache-2.0 WITH LLVM-exception)
// so the harness compiles outside of an LLVM/libfuzzer toolchain — in
// particular under GCC, where the libfuzzer headers aren't shipped.
#include "FuzzedDataProvider.h"

using namespace OpenBabel;

static const std::vector<std::string> &writableFormats() {
    static const std::vector<std::string> formats = []() {
        std::vector<std::string> ids;
        OBPlugin::ListAsVector("formats", "ids", ids);
        std::vector<std::string> writable;
        writable.reserve(ids.size());
        OBConversion probe;
        for (const auto &id : ids) {
            if (probe.SetOutFormat(id.c_str()))
                writable.push_back(id);
        }
        std::sort(writable.begin(), writable.end());
        return writable;
    }();
    return formats;
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    obErrorLog.StopLogging();

    FuzzedDataProvider fdp(Data, Size);

    OBMol mol;
    if (fdp.ConsumeBool())
        mol.SetTitle(fdp.ConsumeRandomLengthString(64).c_str());
    if (fdp.ConsumeBool()) {
        // Some writers emit OBPairData as SDF tags / CML properties / etc.
        OBPairData *pd = new OBPairData;
        pd->SetAttribute(fdp.ConsumeRandomLengthString(16));
        pd->SetValue(fdp.ConsumeRandomLengthString(32));
        mol.SetData(pd);
    }

    for (const auto &fmt : writableFormats()) {
        OBConversion conv;
        if (!conv.SetOutFormat(fmt.c_str()))
            continue;
        try {
            (void)conv.WriteString(&mol, false);
        } catch (...) {
            // writers shouldn't throw, but swallow if they do — we're
            // hunting crashes and leaks, not exception correctness
        }
    }

    return 0;
}

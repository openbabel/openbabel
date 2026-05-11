// Standalone driver for libFuzzer-style targets, used outside OSS-Fuzz.
//
// Each argv entry after argv[0] is treated as a path to a corpus file;
// the file is read into memory and handed to LLVMFuzzerTestOneInput.
// This lets us replay saved OSS-Fuzz reproducers as plain CTest entries
// so the harness doubles as a regression test for previously-fixed bugs.
//
// When LIB_FUZZING_ENGINE is set (i.e. an OSS-Fuzz build), this file is
// not compiled — the fuzzing engine provides its own main.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

static int RunOne(const char *path)
{
  FILE *f = std::fopen(path, "rb");
  if (!f) {
    std::fprintf(stderr, "standalone_main: cannot open %s\n", path);
    return 1;
  }
  std::fseek(f, 0, SEEK_END);
  long sz = std::ftell(f);
  std::fseek(f, 0, SEEK_SET);
  if (sz < 0) {
    std::fclose(f);
    return 1;
  }
  std::vector<std::uint8_t> buf(static_cast<size_t>(sz));
  if (sz > 0 && std::fread(buf.data(), 1, buf.size(), f) != buf.size()) {
    std::fclose(f);
    std::fprintf(stderr, "standalone_main: short read on %s\n", path);
    return 1;
  }
  std::fclose(f);
  LLVMFuzzerTestOneInput(buf.data(), buf.size());
  return 0;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <corpus-file> [<corpus-file>...]\n",
                 argv[0]);
    return 1;
  }
  for (int i = 1; i < argc; ++i) {
    int rc = RunOne(argv[i]);
    if (rc != 0)
      return rc;
  }
  return 0;
}

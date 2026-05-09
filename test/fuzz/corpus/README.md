# Fuzz Corpus

Saved inputs for the libFuzzer harnesses in [`test/fuzz/`](..). Each
subdirectory matches a fuzz-target name; files in it are passed to that
target either by the OSS-Fuzz fuzzer (as seed inputs) or by the
standalone CTest driver (as regression inputs that must not crash).

## Layout

```
corpus/
  fuzz_convert/                 # generic in→out conversion harness
  fuzz_obconversion_sdf/        # ReadString with SetInFormat("sdf")
  fuzz_obconversion_smiles/     # ReadString with SetInFormat("smiles")
```

## File naming

- `clusterfuzz-testcase-minimized-<target>-<id>` — the original OSS-Fuzz
  filename, kept verbatim. The numeric ID is the OSS-Fuzz testcase ID
  and is the canonical way to look up the original report.
- Files added by maintainers as seed corpus may use any short name.

The clusterfuzz files are minimized binary blobs — for `fuzz_convert`
they encode `FuzzedDataProvider` choices (input format, output format,
body), so they are not human-readable chemistry files. To replay one
manually:

```sh
./fuzz_convert path/to/clusterfuzz-testcase-minimized-fuzz_convert-1234
```

## Use as OSS-Fuzz seed corpus

The OSS-Fuzz `build.sh` for openbabel can zip these directories as
seed corpora next to each target binary, e.g.:

```sh
for t in fuzz_convert fuzz_obconversion_sdf fuzz_obconversion_smiles; do
  zip -j "$OUT/${t}_seed_corpus.zip" "$SRC/openbabel/test/fuzz/corpus/$t"/*
done
```

## Adding a new corpus file

1. Drop the file in the matching `corpus/<target>/` directory.
2. Re-run CMake (the file glob is evaluated at configure time).
3. The new file is automatically picked up by the
   `fuzz_<target>_<filename>` regression test.

For raw chemistry-format reproducers (CIF / MOL2 / CDXML / etc. that
trigger a parser bug without going through `FuzzedDataProvider`), see
[`test/files/fuzz_regress/`](../../files/fuzz_regress/) and the
`fuzzregress` test driver instead.

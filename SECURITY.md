# Security Policy

Open Babel is a cheminformatics library, primarily used to read and write
chemistry file formats. The vast majority of inputs come from trusted
sources (a researcher's own files, output from a calculation they ran,
reference databases). Even so, Open Babel is shipped by Linux
distributions and embedded in services that may parse untrusted input,
so we treat memory-safety bugs in format parsers as security issues.

## Reporting a Vulnerability

Please report security vulnerabilities **privately** via GitHub's
"Report a vulnerability" button on the
[Security tab](https://github.com/openbabel/openbabel/security/advisories/new).
This opens a private draft advisory only visible to maintainers.

If you cannot use GitHub, email the maintainers via the address listed
on the [Open Babel website](https://openbabel.org/) and include
`SECURITY` in the subject line. Please do **not** open a public issue
for an unfixed vulnerability.

When reporting, please include:

- A minimized reproducer (input file) and the command used to trigger
  the crash, e.g. `obabel -icif repro.cif -osmi`
- The Open Babel version and commit hash
- Build configuration (compiler, sanitizers, OS)
- Sanitizer output if available (ASAN/UBSAN/MSAN report)

We aim to acknowledge reports within one week and to ship a fix in the
next minor release. Coordinated disclosure timelines are negotiated
case-by-case.

## Supported Versions

Security fixes are applied to the `master` branch and included in the
next tagged release. We do not generally maintain long-term support
branches; downstream packagers (Debian, Fedora, conda-forge, etc.)
are responsible for backports to their own stable releases.

| Version | Supported          |
| ------- | ------------------ |
| 3.1.x   | :white_check_mark: |
| < 3.1   | :x:                |

## CVE Tracking

The table below tracks publicly-assigned CVEs for Open Babel. Items
without a "Fixed in" version remain open. The patch column points to
the commit that lands the fix on `master`; the GHSA column links to
the GitHub Security Advisory once published.

### 2026

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2026-2704](https://nvd.nist.gov/vuln/detail/CVE-2026-2704) | CIF: `transform3d::DescribeAsString` | Out-of-bounds read | Pending merge (#2862) | — | — | — |
| [CVE-2026-2705](https://nvd.nist.gov/vuln/detail/CVE-2026-2705) | MOL2: `OBAtom::SetFormalCharge` | NULL dereference | Pending merge (#2862) | — | — | — |
| [CVE-2026-3408](https://nvd.nist.gov/vuln/detail/CVE-2026-3408) | CDXML: `OBAtom::GetExplicitValence` | NULL dereference | Pending merge (#2862) | — | — | — |

### 2025

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2025-10994](https://nvd.nist.gov/vuln/detail/CVE-2025-10994) | GAMESS: `GAMESSOutputFormat::ReadMolecule` | Use-after-free | Pending merge ([#2834](https://github.com/openbabel/openbabel/issues/2834)) | — | — | — |
| [CVE-2025-10995](https://nvd.nist.gov/vuln/detail/CVE-2025-10995) | zipstream: `basic_unzip_streambuf::underflow` | Overlapping memcpy | Pending merge ([#2832](https://github.com/openbabel/openbabel/issues/2832)) | — | — | — |
| [CVE-2025-10996](https://nvd.nist.gov/vuln/detail/CVE-2025-10996) | SMILES: `OBSmilesParser::ParseSmiles` | Heap-buffer-overflow | Fixed (b34cd604) ([#2831](https://github.com/openbabel/openbabel/issues/2831)) | — | b34cd604 | — |
| [CVE-2025-10997](https://nvd.nist.gov/vuln/detail/CVE-2025-10997) | ChemKin: `ChemKinFormat::CheckSpecies` | Heap-buffer-overflow | Open ([#2830](https://github.com/openbabel/openbabel/issues/2830)) | — | — | — |
| [CVE-2025-10998](https://nvd.nist.gov/vuln/detail/CVE-2025-10998) | ChemKin: `ChemKinFormat::ReadReactionQualifierLines` | NULL dereference | Open ([#2829](https://github.com/openbabel/openbabel/issues/2829)) | — | — | — |
| [CVE-2025-10999](https://nvd.nist.gov/vuln/detail/CVE-2025-10999) | CACAO: `CacaoFormat::SetHilderbrandt` | NULL dereference | Fixed (ecaed96f) ([#2827](https://github.com/openbabel/openbabel/issues/2827)) | — | ecaed96f | — |
| [CVE-2025-11000](https://nvd.nist.gov/vuln/detail/CVE-2025-11000) | PQS: `lowerit` pre-buffer read | Out-of-bounds read | Fixed (duplicate of OSS-Fuzz `lowerit` fix, `f4a5ebae`) | — | f4a5ebae | — |

### 2022 (Cisco TALOS batch — see [#2650](https://github.com/openbabel/openbabel/issues/2650))

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2022-37331](https://nvd.nist.gov/vuln/detail/CVE-2022-37331) | Gaussian: `coords_type` orientation | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-41793](https://nvd.nist.gov/vuln/detail/CVE-2022-41793) | CSR: `PadString` title | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-42885](https://nvd.nist.gov/vuln/detail/CVE-2022-42885) | GRO: res | Uninitialized pointer | Pending merge | — | — | — |
| [CVE-2022-43467](https://nvd.nist.gov/vuln/detail/CVE-2022-43467) | PQS: coord_file | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-43607](https://nvd.nist.gov/vuln/detail/CVE-2022-43607) | MOL2: attribute/value | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-44451](https://nvd.nist.gov/vuln/detail/CVE-2022-44451) | MSI: atom | Uninitialized pointer | Pending merge | — | — | — |
| [CVE-2022-46280](https://nvd.nist.gov/vuln/detail/CVE-2022-46280) | PQS: pFormat | Uninitialized pointer | Pending merge | — | — | — |
| [CVE-2022-46289](https://nvd.nist.gov/vuln/detail/CVE-2022-46289) | ORCA: nAtoms | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46290](https://nvd.nist.gov/vuln/detail/CVE-2022-46290) | ORCA: nAtoms | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46291](https://nvd.nist.gov/vuln/detail/CVE-2022-46291) | translationVectors (Gaussian) | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46292](https://nvd.nist.gov/vuln/detail/CVE-2022-46292) | translationVectors (MOPAC: UNIT CELL TRANSLATION) | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46293](https://nvd.nist.gov/vuln/detail/CVE-2022-46293) | translationVectors (MOPAC: FINAL POINT) | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46294](https://nvd.nist.gov/vuln/detail/CVE-2022-46294) | translationVectors (MOPAC IN: Tv atom) | Out-of-bounds write | Pending merge | — | — | — |
| [CVE-2022-46295](https://nvd.nist.gov/vuln/detail/CVE-2022-46295) | translationVectors (MSI) | Out-of-bounds write | Pending merge | — | — | — |

## Regression Tests

Once a CVE is fixed, the minimized reproducer is checked in under
[`test/files/fuzz_regress/`](test/files/fuzz_regress/) using the
naming convention `cve-YYYY-NNNNN.<ext>`. The
[`fuzzregresstest`](test/fuzzregresstest.cpp) harness re-runs each
reproducer through `OBConversion::ReadFile` on every CI build so that
the original crash cannot return undetected. CI builds at least one
job with `-fsanitize=address,undefined` so memory-safety regressions
are caught immediately.

## Threat Model and Scope

- **In scope:** memory-safety bugs (heap/stack overflows, use-after-free,
  uninitialized reads, NULL dereferences) reachable through the public
  `OBConversion::ReadFile` / `WriteFile` API, the `obabel` command-line
  tool, or the language bindings.
- **Out of scope by default:** denial-of-service from very large or
  pathological inputs (we accept that parsing untrusted input may be
  slow); incorrect chemistry on malformed inputs that does not involve
  memory unsafety; bugs in third-party libraries that ship as
  dependencies (please report those upstream as well).

If you are unsure whether something qualifies, report it anyway and
let us decide.

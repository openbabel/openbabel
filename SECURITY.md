# Security Policy

Open Babel is a cheminformatics library, primarily used to read and write
chemistry file formats. The vast majority of inputs come from trusted
sources (a researcher's own files, output from a calculation they ran,
reference databases) for local use. Even so, Open Babel is shipped by Linux
distributions and embedded in services that may parse untrusted input,
so we treat memory-safety bugs in format parsers as potential security
issues and evaluate their impact based on *exploitability* and realistic
usage scenarios.

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

We aim to acknowledge reports within one to two weeks and to ship a fix 
in the next minor release. Coordinated disclosure timelines are negotiated
case-by-case.

Please note that not all reported issues will be classified as security
vulnerabilities or assigned CVEs; classification is based on impact and
exploitability.

## Supported Versions

Security fixes are applied to the `master` branch and included in the
next tagged release. We do not generally maintain long-term support
branches; downstream packagers (Debian, Fedora, conda-forge, etc.)
are responsible for backports to their own stable releases.

| Version | Supported          |
| ------- | ------------------ |
| 3.2.x   | :white_check_mark: |
| < 3.2   | :x:                |

## CVE Tracking

The table below tracks publicly-assigned CVEs for Open Babel. Items
without a "Fixed in" version remain open. The patch column points to
the commit that lands the fix on `master`; the GHSA column links to
the GitHub Security Advisory once published.

### 2026

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2026-54751] | MMD `MacroModFormat::ReadMolecule` | Uninitialized read | Fixed (a08065d) ([#2973](https://github.com/openbabel/openbabel/pull/2973)) | 3.2.1 | a08065d | [GHSA-gxqq-gwf6-3xmq](https://github.com/openbabel/openbabel/security/advisories/GHSA-gxqq-gwf6-3xmq) |
| [CVE-2026-2704](https://nvd.nist.gov/vuln/detail/CVE-2026-2704) | CIF: `transform3d::DescribeAsString` | Out-of-bounds read | Fixed (e23a224b) ([#2862](https://github.com/openbabel/openbabel/pull/2862)) | 3.2.0 | e23a224b | [GHSA-6xw4-2g22-26h8](https://github.com/openbabel/openbabel/security/advisories/GHSA-6xw4-2g22-26h8) |
| [CVE-2026-2705](https://nvd.nist.gov/vuln/detail/CVE-2026-2705) | MOL2: `OBAtom::SetFormalCharge` | NULL dereference | Fixed (e23a224b) ([#2862](https://github.com/openbabel/openbabel/pull/2862)) | 3.2.0 | e23a224b | [GHSA-4w5w-4fhm-q483](https://github.com/openbabel/openbabel/security/advisories/GHSA-4w5w-4fhm-q483) |
| [CVE-2026-3408](https://nvd.nist.gov/vuln/detail/CVE-2026-3408) | CDXML: `OBAtom::GetExplicitValence` | NULL dereference | Fixed (e23a224b) ([#2862](https://github.com/openbabel/openbabel/pull/2862)) | 3.2.0 | e23a224b | [GHSA-rxpr-wq63-jr7p](https://github.com/openbabel/openbabel/security/advisories/GHSA-rxpr-wq63-jr7p) |

### 2025

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2025-10994](https://nvd.nist.gov/vuln/detail/CVE-2025-10994) | GAMESS: `GAMESSOutputFormat::ReadMolecule` | Use-after-free | Fixed (95033d27) ([#2834](https://github.com/openbabel/openbabel/issues/2834)) | 3.2.0 | 95033d27 | [GHSA-pp85-5j63-xpq3](https://github.com/openbabel/openbabel/security/advisories/GHSA-pp85-5j63-xpq3) |
| [CVE-2025-10995](https://nvd.nist.gov/vuln/detail/CVE-2025-10995) | zipstream: `basic_unzip_streambuf::underflow` | Overlapping memcpy | Fixed (d4621d41) ([#2832](https://github.com/openbabel/openbabel/issues/2832)) | 3.2.0 | d4621d41 | [GHSA-8j3x-m868-cpw8](https://github.com/openbabel/openbabel/security/advisories/GHSA-8j3x-m868-cpw8) |
| [CVE-2025-10996](https://nvd.nist.gov/vuln/detail/CVE-2025-10996) | SMILES: `OBSmilesParser::ParseSmiles` | Heap-buffer-overflow | Fixed (b34cd604) ([#2831](https://github.com/openbabel/openbabel/issues/2831)) | 3.2.0 | b34cd604 | [GHSA-j35x-w4gj-pf7w](https://github.com/openbabel/openbabel/security/advisories/GHSA-j35x-w4gj-pf7w) |
| [CVE-2025-10997](https://nvd.nist.gov/vuln/detail/CVE-2025-10997) | ChemKin: `ChemKinFormat::CheckSpecies` | Heap-buffer-overflow | Fixed (af4a4212) ([#2830](https://github.com/openbabel/openbabel/issues/2830)) | 3.2.0 | af4a4212 | [GHSA-8wq6-qh76-wpv9](https://github.com/openbabel/openbabel/security/advisories/GHSA-8wq6-qh76-wpv9) |
| [CVE-2025-10998](https://nvd.nist.gov/vuln/detail/CVE-2025-10998) | ChemKin: `ChemKinFormat::ReadReactionQualifierLines` | NULL dereference | Fixed (af4a4212) ([#2829](https://github.com/openbabel/openbabel/issues/2829)) | 3.2.0 | af4a4212 | [GHSA-j9pp-7wfg-q7fj](https://github.com/openbabel/openbabel/security/advisories/GHSA-j9pp-7wfg-q7fj) |
| [CVE-2025-10999](https://nvd.nist.gov/vuln/detail/CVE-2025-10999) | CACAO: `CacaoFormat::SetHilderbrandt` | NULL dereference | Fixed (ecaed96f) ([#2827](https://github.com/openbabel/openbabel/issues/2827)) | 3.2.0 | ecaed96f | [GHSA-55j6-rjhx-hwfh](https://github.com/openbabel/openbabel/security/advisories/GHSA-55j6-rjhx-hwfh) |
| [CVE-2025-11000](https://nvd.nist.gov/vuln/detail/CVE-2025-11000) | PQS: `lowerit` pre-buffer read | Out-of-bounds read | Fixed (duplicate of OSS-Fuzz `lowerit` fix, `f4a5ebae`) | 3.2.0 | f4a5ebae | [GHSA-m982-7q3h-r784](https://github.com/openbabel/openbabel/security/advisories/GHSA-m982-7q3h-r784) |

### 2022 (Cisco TALOS batch — see [#2650](https://github.com/openbabel/openbabel/issues/2650))

| CVE | Component | Type | Status | Fixed in | Patch | GHSA |
| --- | --- | --- | --- | --- | --- | --- |
| [CVE-2022-37331](https://nvd.nist.gov/vuln/detail/CVE-2022-37331) | Gaussian: `coords_type` orientation | Out-of-bounds write | Fixed (528c142f) | 3.2.0 | 528c142f | [GHSA-vr3p-gg26-45v9](https://github.com/openbabel/openbabel/security/advisories/GHSA-vr3p-gg26-45v9) |
| [CVE-2022-41793](https://nvd.nist.gov/vuln/detail/CVE-2022-41793) | CSR: `PadString` title | Out-of-bounds write | Fixed (528c142f) | 3.2.0 | 528c142f | [GHSA-p594-7xw4-g76p](https://github.com/openbabel/openbabel/security/advisories/GHSA-p594-7xw4-g76p) |
| [CVE-2022-42885](https://nvd.nist.gov/vuln/detail/CVE-2022-42885) | GRO: res | Uninitialized pointer | Fixed (fa9a2d9a) | 3.2.0 | fa9a2d9a | [GHSA-mw5r-wq2m-397c](https://github.com/openbabel/openbabel/security/advisories/GHSA-mw5r-wq2m-397c) |
| [CVE-2022-43467](https://nvd.nist.gov/vuln/detail/CVE-2022-43467) | PQS: coord_file | Out-of-bounds write | Fixed (2a7d2cda) | 3.2.0 | 2a7d2cda | [GHSA-f29h-2h58-48r7](https://github.com/openbabel/openbabel/security/advisories/GHSA-f29h-2h58-48r7) |
| [CVE-2022-43607](https://nvd.nist.gov/vuln/detail/CVE-2022-43607) | MOL2: attribute/value | Out-of-bounds write | Fixed (4110d59a) | 3.2.0 | 4110d59a | [GHSA-vjg6-gm8m-v5g6](https://github.com/openbabel/openbabel/security/advisories/GHSA-vjg6-gm8m-v5g6) |
| [CVE-2022-44451](https://nvd.nist.gov/vuln/detail/CVE-2022-44451) | MSI: atom | Uninitialized pointer | Fixed (fa9a2d9a) | 3.2.0 | fa9a2d9a | [GHSA-jr2x-6qf6-q5mc](https://github.com/openbabel/openbabel/security/advisories/GHSA-jr2x-6qf6-q5mc) |
| [CVE-2022-46280](https://nvd.nist.gov/vuln/detail/CVE-2022-46280) | PQS: pFormat | Uninitialized pointer | Fixed (2a7d2cda) | 3.2.0 | 2a7d2cda | [GHSA-8qxc-57hf-hc9j](https://github.com/openbabel/openbabel/security/advisories/GHSA-8qxc-57hf-hc9j) |
| [CVE-2022-46289](https://nvd.nist.gov/vuln/detail/CVE-2022-46289) | ORCA: nAtoms | Out-of-bounds write | Fixed (b239d06e) | 3.2.0 | b239d06e | [GHSA-rj4c-r689-cm87](https://github.com/openbabel/openbabel/security/advisories/GHSA-rj4c-r689-cm87) |
| [CVE-2022-46290](https://nvd.nist.gov/vuln/detail/CVE-2022-46290) | ORCA: nAtoms | Out-of-bounds write | Fixed (b239d06e) | 3.2.0 | b239d06e | [GHSA-5rff-8f7c-8jmw](https://github.com/openbabel/openbabel/security/advisories/GHSA-5rff-8f7c-8jmw) |
| [CVE-2022-46291](https://nvd.nist.gov/vuln/detail/CVE-2022-46291) | translationVectors (Gaussian) | Out-of-bounds write | Fixed (40e85213) | 3.2.0 | 40e85213 | [GHSA-jg3h-pv7c-4f9c](https://github.com/openbabel/openbabel/security/advisories/GHSA-jg3h-pv7c-4f9c) |
| [CVE-2022-46292](https://nvd.nist.gov/vuln/detail/CVE-2022-46292) | translationVectors (MOPAC: UNIT CELL TRANSLATION) | Out-of-bounds write | Fixed (40e85213) | 3.2.0 | 40e85213 | [GHSA-55f6-pf8r-c2f4](https://github.com/openbabel/openbabel/security/advisories/GHSA-55f6-pf8r-c2f4) |
| [CVE-2022-46293](https://nvd.nist.gov/vuln/detail/CVE-2022-46293) | translationVectors (MOPAC: FINAL POINT) | Out-of-bounds write | Fixed (40e85213) | 3.2.0 | 40e85213 | [GHSA-7h6r-6p76-68c9](https://github.com/openbabel/openbabel/security/advisories/GHSA-7h6r-6p76-68c9) |
| [CVE-2022-46294](https://nvd.nist.gov/vuln/detail/CVE-2022-46294) | translationVectors (MOPAC IN: Tv atom) | Out-of-bounds write | Fixed (40e85213) | 3.2.0 | 40e85213 | [GHSA-mjmg-352j-f456](https://github.com/openbabel/openbabel/security/advisories/GHSA-mjmg-352j-f456) |
| [CVE-2022-46295](https://nvd.nist.gov/vuln/detail/CVE-2022-46295) | translationVectors (MSI) | Out-of-bounds write | Fixed (40e85213) | 3.2.0 | 40e85213 | [GHSA-f8h2-c479-vqxf](https://github.com/openbabel/openbabel/security/advisories/GHSA-f8h2-c479-vqxf) |

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

- **In scope:** High-impact memory-safety bugs (e.g., Remote Code Execution vectors, arbitrary heap/stack overflows, use-after-free) reachable through the public `OBConversion::ReadFile` / `WriteFile` API under realistic invocation.

- **Out of scope by default:** 
  - **Standard application crashes**, assertions, or unhandled exceptions resulting in a simple Denial of Service (DoS) for a local CLI run.
  - **Memory leaks** or minor uninitialized reads that do not impact execution flow.
  - **Individual NULL pointer dereferences** in legacy or rarely used file formats; these will be treated as standard application bugs rather than security vulnerabilities unless a clear exploit vector is demonstrated.
  - *Denial-of-service** from very large or
  pathological inputs (we accept that parsing untrusted input may be
  slow); incorrect chemistry on malformed inputs that does not involve
  memory unsafety; bugs in third-party libraries that ship as
  dependencies (please report those upstream as well).

Denial-of-service issues (e.g., excessive CPU or memory usage) are
generally not considered security vulnerabilities unless they provably affect
long-running services or network-exposed deployments.

If you are unsure whether something qualifies, report it anyway and
let us decide.

## Impact Classification

Not all bugs reachable via file parsing are treated as security
vulnerabilities.

In general:

- **Security vulnerabilities** include memory corruption issues with a
  plausible path to code execution, data corruption, or other impact
  beyond a crash, especially when parsing untrusted input.

- **Bug fixes (no CVE)** include:
  - NULL dereferences leading only to crashes
  - Out-of-bounds reads without control over program flow
  - Errors requiring highly contrived or non-realistic inputs
  - Issues unlikely to occur in real-world usage
  - Bugs in bundled third-party code (please report upstream)

Open Babel is primarily a local-use library and CLI tool. Many issues
require a user to explicitly process a malicious file; absent a realistic
exploitation scenario, these are typically treated as correctness bugs
rather than security vulnerabilities.

We may still fix and credit such reports, but they do not generally
receive CVE assignments.

If you believe a crash-class bug is actually a memory-corruption path, escalate it with a reproducer and an ASAN/UBSAN trace showing the out-of-bounds write or use-after-free. We take those seriously regardless of how unlikely exploitation seems.

## Bulk Reports and CVE Consolidation Policy

To prevent alert fatigue for downstream package maintainers, Open Babel enforces a strict consolidation policy for security advisories:
1. **One Advisory per Batch:** If an automated tool or fuzzer discovers multiple edge-case parsing bugs (e.g., multiple NULL dereferences across different formats), they **must** be submitted as a single consolidated report, or they will be merged by the maintainers.
2. **Consolidated CVEs:** Open Babel will issue at most **one CVE per minor release batch** for clusters of low-severity parser bugs. We will not issue individual CVEs for multiple distinct files fixed in the same patch lifecycle.
3. **Credit:** Reporters will still receive full credit in the release notes and the consolidated GHSA/CVE for all discovered bugs.
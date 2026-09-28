# Security review - 2026-09-28

Scope: ITensorSchmidtPlots 0.2.0 source, examples, tests, Git history, GitHub
workflow, and the resolved Julia 1.11.5 dependency environment on macOS arm64.
This is a source review and advisory comparison, not a penetration test or a
guarantee that the entire Julia and native-library stack is vulnerability-free.

## Assessment and changes

No exploitable vulnerability was identified in this package's own source during
the review. Calculations use in-memory numeric inputs; PDF export writes to a
caller-selected path. There is no network request, shell execution, dynamic
evaluation, deserialization, external LaTeX compiler, or input-dependent source
include in the package. Malformed numeric input is rejected, and the MPS is
copied before gauge changes. Arbitrary Julia objects and custom Makie functions
are trusted code, not sandboxed input.

Dependency updates removed **27 advisory matches** from the plotting libraries:

| Library | Previously resolved | Updated and tested | Matches removed |
| --- | --- | --- | ---: |
| FFMPEG_jll | 8.1.2+0 | 9.0.1+0 | 22 |
| Giflib_jll | 5.2.3+0 | 6.1.3+0 | 3 |
| HarfBuzz_jll | 8.5.1+0 | 100.14004.0+0 | 2 |

Project.toml now requires CairoMakie >= 0.15.15 and Makie >= 0.24.15 within their
compatible minor releases. It also explicitly requires FFMPEG_jll >= 9.0.1,
Giflib_jll >= 6.1.3, and HarfBuzz_jll >= 100.14002.1 within their compatible major
releases. These already-transitive libraries are declared directly so that
users' existing environments cannot silently keep the affected older versions.
Related Pango, libass, libwebp, and other compatible dependencies were updated.

The full [resolved version inventory](security-dependencies.csv) includes
283 packages, 282 of which carry explicit versions in the manifest. This snapshot
is specific to the reviewed environment; Julia may resolve different versions
on other platforms or in existing projects.

## Advisory comparison and remaining findings

The scan used all **1,717 published advisories** in the official
[Julia advisory database at dcb1dba](https://github.com/JuliaLang/SecurityAdvisories.jl/tree/dcb1dba21affb432702e3ec0e08d5753343c5a10).
That was the repository's current revision when fetched; its commit date is
2026-08-24. This database cannot establish coverage of every later disclosure.
Withdrawn entries were excluded. Package names and affected ranges were compared
against the manifest using Julia's version ordering, including JLL build suffixes.
Unrecognized range syntax causes the scanner to stop, rather than skip an entry.
The database explicitly notes that native-artifact coverage is incomplete.

The initial environment produced 106 advisory matches; the updated environment
produces **79**. These are version-range matches, not 79 demonstrated exploits
of this package. The [complete remaining matches](security-advisory-matches.csv)
are retained for inspection rather than suppressed.

| Remaining component | Version | Matches | Assessment |
| --- | --- | ---: | --- |
| LibCURL_jll | 8.6.0+0 | 37 | Bundled with the old Julia runtime; networking advisories. |
| MbedTLS_jll | 2.28.6+0 | 17 | Bundled runtime TLS library. |
| LibSSH2_jll | 1.11.0+1 | 11 | Bundled runtime SSH library. |
| LibGit2_jll | 1.7.2+0 | 6 | Bundled runtime Git library. |
| nghttp2_jll | 1.59.0+0 | 3 | Bundled runtime HTTP/2 library; some entries concern the proxy executable. |
| Zlib_jll | 1.2.13+1 | 2 | Bundled runtime compression library; indirect use is possible. |
| p7zip_jll | 17.4.0+2 | 1 | Advisory concerns Windows archive extraction; not demonstrated on this macOS PDF path. |
| Cairo_jll | 1.18.7+0 | 1 | JLSEC-2026-1194 describes PDF-content memory remanence; no fixed range is published. |
| Libmount_jll | 2.42.0+0 | 1 | JLSEC-2026-1263 concerns libblkid partition probing; this package does not probe or mount disks. |

**Action for users:** update Julia itself to a currently maintained stable or LTS
release, then update the project dependencies and repeat the scan. `Pkg.update()`
alone does not replace these bundled runtime libraries. The current stable release
was [Julia 1.13.1](https://julialang.org/downloads/manual-downloads/) on the review
date; the local regression tests used the existing Julia 1.11.5 installation.
No claim is made that merely upgrading Julia resolves every remaining advisory.

The package has no networking or media-decoding feature. Its presets use bundled
fonts and generated vector geometry. Those properties reduce the relevant attack
surface, but do not establish complete transitive call-graph reachability proofs.
Do not treat Julia or this package as a sandbox, or assume sensitive arrays and
PDF contents are securely erased from process memory. See [SECURITY.md](../SECURITY.md).

## Repository and workflow checks

- Scanned 30 pre-release history blobs and current source for common private-key,
  GitHub-token, AWS-key, and API-secret patterns; no matches found. Pattern scans
  cannot identify every possible secret.
- Reviewed CI: read-only `contents` permission, immutable action commit pins,
  no persisted checkout credentials, no secrets, no `pull_request_target`,
  no privileged/self-hosted runner, and a 30-minute job timeout.
- The public advisory APIs for `actions/checkout` and `julia-actions/setup-julia`
  returned no published repository advisories on the review date. This does not
  audit their entire dependency trees.
- Checked the current public advisory APIs for `ITensor/ITensors.jl`,
  `ITensor/ITensorMPS.jl`, `MakieOrg/Makie.jl` (including CairoMakie), and
  `JuliaStrings/LaTeXStrings.jl`; all returned empty lists. Their repository
  identities were obtained from Julia's refreshed General registry.
- All 846 regression tests passed with the updated dependencies. They cover
  numerical correctness, invalid inputs, preservation of
  the supplied MPS, truncation boundaries, native Makie composition, and PDF
  output. The removed preset is rejected before writing an output file.
- Rebuilt all three example PDFs, verified 156 truncation ranks independently
  against dense SVDs, and inspected the rendered figures and embedded fonts.

## Reproduce the dependency report

Instantiate the project to produce a Manifest.toml, obtain the official advisory
database, and run the included report-only scanner:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
git clone https://github.com/JuliaLang/SecurityAdvisories.jl.git /path/to/advisories
julia --startup-file=no security/audit_dependencies.jl Manifest.toml /path/to/advisories report.toml
```

For the historical comparison, check out database commit
`dcb1dba21affb432702e3ec0e08d5753343c5a10`. For a current review, use current
published advisories and the manifest of the environment actually used.
The scanner does not execute advisory content or send the manifest anywhere.
It reports matching ranges and exits normally when report generation succeeds;
that exit status is **not** a clean bill of health or a CI security gate.

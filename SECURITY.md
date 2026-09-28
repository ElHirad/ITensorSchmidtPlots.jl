# Security policy

Use the latest package release and the latest patch of a maintained Julia stable
or LTS release. Julia version compatibility alone does not imply that an old
runtime still receives security fixes.

## Intended use

This package computes spectra from in-memory MPSs or numeric vectors and writes
PDFs to a caller-selected path. It provides no network service, remote-data
loader, shell command interface, or external LaTeX execution. The bundled
Computer Modern fonts are used by the preset and examples.

Julia code and custom Makie functions run with the caller's permissions. This
library is not an isolation boundary for untrusted Julia objects, executable
input, fonts, or uploads. Large MPSs can consume substantial memory and CPU. If
you build a service around it, enforce resource and output-path restrictions in
that service, in a separate process with appropriate permissions.

`save_schmidt` accepts only PDF filenames. It replaces an existing file at the
requested path, following Julia/Makie filesystem semantics, including symlinks.
Supply a path you control. Spectrum extraction copies the input MPS.

## Dependencies

Version 0.2.0 enforces patched minimum versions of the already-transitive
FFMPEG_jll, Giflib_jll, and HarfBuzz_jll libraries. Julia's resolver must satisfy
these bounds; do not bypass them with locally modified dependencies.

The [dated review](docs/SECURITY_REVIEW.md) records the tested versions, security
checks, fixes, remaining advisory matches, and reproducible scan instructions.
An absence of findings is not a guarantee against unknown vulnerabilities.

## Reporting

Report potentially exploitable issues privately to the maintainer or the
affected upstream project according to its security policy. If no private
maintainer contact is available, the [Julia security team](https://security.julialang.org/about/)
accepts Julia ecosystem reports at `security@julialang.org`.
Do not put credentials, private datasets, or exploit details in a public issue.
Include the package and Julia versions, a minimal reproduction, the dependency
manifest, and the affected operation; remove private paths and data first.

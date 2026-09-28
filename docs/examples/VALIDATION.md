# Numerical validation

Generated with Julia 1.11.5, ITensors 0.9.30, ITensorMPS 0.3.45, CairoMakie 0.15.13.
Every internal bond is compared with a direct dense SVD of the original vector.
Truncation ranks are independently checked against dense-SVD squared tails.

| Example | Max absolute SVD error | Analytic error | Largest tail error | Rank checks |
| --- | ---: | ---: | ---: | ---: |
| Nonuniform entangled pairs | 8.326672684688674e-16 | 8.881784197001252e-16 | 3.9898639947466563e-17 | 52 / 52 |
| Synthetic multiscale field | 1.5543122344752192e-15 | — | 9.194034422677078e-17 | 52 / 52 |
| Complex random state | 2.4424906541753444e-15 | — | 2.6020852139652106e-17 | 52 / 52 |

The first example additionally agrees with the analytic product of the seven pair spectra.
Values are normalized singular values (not squared probabilities). Tests also cover complex amplitudes, QN blocks, gauge changes, input preservation, zero states, rank caps, and PDF export.
These are synthetic examples, not reconstructions of the papers' DNS data.

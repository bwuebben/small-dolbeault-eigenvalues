# Computational checks

Each script uses only the Python 3 standard library and exits with an
error if any assertion fails. Run one from the repository root, for
example

```bash
python3 checks/check_two_step_flow.py
```

or all seven with `make checks`. Exact identities use rational
arithmetic (`fractions.Fraction`); spectra are computed by Jacobi
diagonalization in 110-digit `decimal` arithmetic. The scripts check the
finite-dimensional algebra and model computations on which the proofs
rely. They do not discretize a bundle or a curve, and they do not verify
the elliptic and parabolic estimates, which are proved in the paper.

| Script | Paper | What it checks |
|---|---|---|
| `check_two_factor_algebra.py` | Section 2 | Block signs, curvature and energy identities in 48 rational cases (ranks 2–8, one to three form components); the sign of the fourth-order coefficient and the sixth-order remainder in a finite-dimensional operator model. |
| `check_several_factor.py` | Sections 3–5 | Exact rational relative-form and Schur-complement bounds with coupled spectral comparisons; 36 three-factor spectra with unequal masses and exponents; the five-factor coefficients 5/6, 3/2, 2, 2; the example showing that the harmonic cancellation cannot be dropped. |
| `check_repeated_factors.py` | Section 13 | The commutator form on multiplicity matrices; the A ⊕ B ⊕ A ⊕ C example with five positive small modes and leading coefficients 1, 1, 1, 2, 2; a nonsimple connected self-extension; interference between entries sharing one harmonic form; 12 coupled relative-form, Schur and kernel checks. |
| `check_harmonic_normal_form.py` | Section 6 | The triangular gauge induction of Lemma 6.1 over Laurent polynomials: preservation of earlier diagonals, the full gauge identity and weighted scaling in 36 exact cases with 2–7 blocks, and a three-block example in which a harmonic (1,3) entry arises from nonharmonic entries. |
| `check_yang_mills_flow.py` | Section 7 | Metric and endomorphism-form comparisons (the Rayleigh exponent four is necessary); the two-vertex reduced flow λ′ = −2λ² with τλ → 1/2 in 100 rank, volume and norm cases; the cubic-residual tail integral used in the late-time comparison. |
| `check_several_factor_flow.py` | Sections 9 and 11 | Moment map, time normalization and energy dissipation in 36 exact cases; the three- and four-factor changes of variables; separated graph spectra; closed scalar solutions from unequal initial data, with edge weights a, c ∼ 1/(4τ) and b ∼ 1/(2τ log τ). |
| `check_two_step_flow.py` | Sections 9 and 11 | Rank- and volume-weighted moment and energy identities; the unequal-rank changes of variables; matrix scaling compared with the log-time model flow; both zigzag balance obstructions; cycle and unequal-mass graph spectra; the coefficient 2/(27 ρ^(2/3)) of Proposition 11.3 at four initial scales. |

`check_several_factor.py` and `check_two_factor_algebra.py` also provide
small matrix helpers imported by the other scripts, so the scripts should
be run from their own directory or by path as above.

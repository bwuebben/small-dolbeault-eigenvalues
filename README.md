# Small Dolbeault eigenvalues along Hermitian Yang–Mills flow

Research manuscript and computational checks by **Bernd Johannes Wuebben**.

**[Read the paper (PDF, 45 pages)](paper/small-dolbeault-eigenvalues.pdf)**, initial version 20 September 2026, this version 30 September 2026.

## Summary

The summand projections of a polystable bundle are holomorphic
endomorphisms, hence zero modes of the Dolbeault Laplacian on
endomorphisms. Along Hermitian Yang–Mills flow a semistable bundle
approaches its graded limit, and the corresponding eigenvalues tend to
zero. The paper determines how fast they vanish, with which leading
coefficients, and which of these coefficients remember the metric the
flow started from. The results cover simple bundles on a curve whose
Jordan–Hölder factors are pairwise nonisomorphic, and extensions of one
polystable bundle by another on compact Kähler manifolds of any dimension.

Throughout, λ₁ ≤ λ₂ ≤ … denote the eigenvalues of the Dolbeault
Laplacian on trace-free endomorphisms, and the flow is normalized as
h⁻¹∂_τh = −2(K_h − κ Id).

## Main results

**Comparison with a finite-dimensional flow (Theorem 9.5).** Put the
bundle in a harmonic normal form adapted to its socle filtration. The
nonzero extension entries define a directed graph with vertex masses
given by the ranks. For every smooth initial metric and every
sufficiently late time T, a solution g^[T] of an explicit flow on
diagonal metrics satisfies

```math
e^{-\delta_T}\,\nu_j\bigl(g^{[T]}(\tau)\bigr)\;\le\;\lambda_j\bigl(h(\tau)\bigr)\;\le\;e^{\delta_T}\,\nu_j\bigl(g^{[T]}(\tau)\bigr),
\qquad \tau\ge T,\qquad \delta_T\to 0 ,
```

where the ν_j are the eigenvalues of the weighted graph Laplacian. The
relative error tends to zero uniformly on the future interval, also for
eigenvalues that decay at different rates.

**The complete leading small spectrum (Theorems 1.1 and 10.1).** The
exponents are governed by the weight grading of Haiden, Katzarkov,
Kontsevich and Pandit, the solution of a quadratic program on the graph.
Suppose some Lagrange multiplier vector of this program is strictly
positive on every active constraint. Then:

- every small eigenvalue satisfies τ^(a_j) λ_j → c_j > 0, with exponents
  determined by the graph and the ranks through a contraction rule;
- the coefficients of the eigenvalues of order 1/τ do not depend on the
  initial metric, and 1/2 is always one of them, with the weight grading
  as eigenvector;
- at each faster rate the coefficients depend on the initial metric, and
  their sum is unbounded over metrics on one fixed bundle;
- the spectral subspaces converge smoothly to explicit eigenspaces of
  weighted graph Laplacians, which identify the stable factors.

**Convergence of the remainder (Theorem 8.3, Corollary 10.2).** Under the
same condition, the bounded remainder in the asymptotics of the diagonal
flow converges, and so does the bounded term in the metric asymptotics
along the actual flow. The limit depends on the initial metric only
through one constant for each connected component of the active
constraints.

**Examples (Theorem 7.1 and Section 11).** For a nonsplit extension of
two stable bundles, τλ₁ → 1/2, independently of the ranks, the extension
class, the initial metric and the dimension. For four stable factors of
ranks (2,1,1,2) in a zigzag, λ₁ ∼ c(h(0)) τ^(−4/3), and explicit initial
metrics h_ρ on one bundle give c(h_ρ) = (2/27) ρ^(−2/3) (1 + O(ρ)). For
four line bundles in the same zigzag, the multiplier condition fails and
λ₁ ∼ 1/(2τ log τ).

**Higher dimensions (Theorem 12.1).** For extensions of one polystable
bundle by another on a compact Kähler manifold, with pairwise
nonisomorphic stable factors, the comparison and the complete-spectrum
theorem hold in every dimension. The input is smooth convergence of the
flow when the graded object is locally free (Sibley and Wentworth).

**Harmonic families with fixed metrics (Sections 2–6, 13–14).** For scaled
harmonic extensions, every small eigenvalue agrees with a weighted graph
Laplacian up to a relative error (Theorem 3.1). A contraction rule
determines all decay orders and limiting subspaces (Theorem 4.1). A
harmonic normal form (Lemma 6.1) constructs such families for every simple
strictly semistable bundle on a curve, including bundles with repeated
stable factors. The limiting endomorphism algebra then recovers the graded
factors and their multiplicities (Theorems 13.1 and 14.1).

**A slope inequality, sharp in every rank (Theorem 1.3, Corollary 15.9).**
For any Hermitian metric h on a bundle of rank r over a compact connected
Kähler manifold, and every saturated subsheaf S of rank s,

```math
\mu(E)-\mu(\mathcal S)\;\ge\;\frac{V(r-s)}{2\pi r}\Bigl(\lambda_1(h)-\frac{r}{\max\{s,\,r-s\}}\,\varepsilon_1(h)\Bigr),
```

where ε₁ is the averaged operator norm of the trace-free mean curvature.
Corrected metrics on nonsplit extensions show that the coefficient
r/max{s, r−s} cannot be decreased for any r and s. The same holds for the
criterion with the supremum norm of the curvature.

## Relation to earlier work

Convergence of the flow on curves is classical (Daskalopoulos; Råde); in
higher dimensions the paper uses Sibley and Wentworth's identification of
the singular set. Haiden, Katzarkov, Kontsevich and Pandit introduced the
weight grading and the weighted-graph flows. Their work contains the
multiplier condition, the exponent 4/3, the logarithmic zigzag, the Green
operator correction and the order comparison, all with a bounded
remainder. Lam extends the metric asymptotics to compact Kähler manifolds,
including repeated factors, also with a bounded remainder. This paper
transfers precise coefficients to the Dolbeault spectrum, proves
convergence of the remainder when the multipliers are positive, and
derives the resulting dichotomy. After a change of variables the diagonal
flow is gradient flow for an exponential loss with separable data. Soudry,
Hoffer, Nacson, Gunasekar and Srebro prove convergence of the remainder
when the support vectors span the data, which corresponds here to a single
connected component of active constraints. Sektnan and Tipler use scaled
harmonic extensions in an adiabatic Hermitian Yang–Mills problem. Precise
references are given in the introduction of the paper.

## Organization

```text
paper/
  small-dolbeault-eigenvalues.pdf  the manuscript
  main.tex                         Section 1 and 16: introduction, open problems
  two_factor.tex                   Section 2: extensions of two stable bundles
  several_factors.tex              Sections 3–5: several summands, decay orders, examples
  normal_form.tex                  Section 6: harmonic normal form on curves
  yang_mills_flow.tex              Section 7: the flow for two factors
  diagonal_flow.tex                Section 8: the diagonal flow and its asymptotics
  transfer.tex                     Section 9: comparison with the diagonal flow
  complete_spectrum.tex            Section 10: the complete leading small spectrum
  flow_examples.tex                Section 11: examples along the flow
  higher_dim.tex                   Section 12: two-step extensions in higher dimensions
  repeated_factors.tex             Section 13: repeated stable summands
  curve_recovery.tex               Section 14: recovery of the graded object
  slope.tex                        Section 15: the slope inequality and its sharpness
checks/                            seven computational checks, indexed in checks/README.md
```

## Build and verify

The checks use only the Python 3 standard library. Building the
manuscript requires TeX Live with `latexmk`. From the repository root:

```bash
make checks   # run the seven computational checks
make paper    # build the manuscript into build/main.pdf
make verify   # run the checks, then build the manuscript
```

The versioned PDF in `paper/` is never overwritten by a build.

The checks verify the finite-dimensional algebra and model computations
used in the proofs: block identities, exact rational quadratic-form
bounds, graph spectra in 110-digit arithmetic, changes of variables in
the reduced flows, and the explicit coefficient 2/(27 ρ^(2/3)). They do
not discretize the bundle and do not verify the elliptic and parabolic
estimates, which are proved in the paper. See the
[guide to the checks](checks/README.md).

## Citation

A [BibTeX entry](CITATION.bib) and machine-readable
[citation metadata](CITATION.cff) are provided.

Code is MIT licensed; the manuscript and repository prose are
CC BY 4.0. See [LICENSE](LICENSE).

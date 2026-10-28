# Small Dolbeault eigenvalues near stability walls

**Bernd Johannes Wuebben** · 27 October 2026 · 13 pages

[Current version (PDF)](paper/main.pdf) · [Manuscript sources](paper/)
· [Abstract and citation](https://bwuebben.github.io/papers/small-dolbeault-stability-walls/)
· [Zenodo:23717764](https://zenodo.org/records/23717764)

Version DOI: [10.5281/zenodo.23717764](https://doi.org/10.5281/zenodo.23717764).
ArXiv submission is pending. [BibTeX](CITATION.bib) and
[citation metadata](CITATION.cff) identify this version.

This paper determines the small endomorphism Dolbeault spectrum of
Hermitian–Einstein bundles as the Kähler class approaches an intersection
of slope-stability walls. It is a separate companion to
[Small Dolbeault eigenvalues along Hermitian Yang–Mills flow](../paper/README.md).
The [repository guide](../README.md) indexes both papers and their separate citations.

## Main results

Let a holomorphic bundle on a compact connected Kähler manifold of dimension
at least two fit into an extension

```math
0\longrightarrow F_0\longrightarrow E\longrightarrow
\bigoplus_{j=1}^kF_j\longrightarrow0.
```

The factors are pairwise nonisomorphic and stable of equal slope at the
limiting class, and every extension component is nonzero. Write
$r_i=\mathrm{rk}F_i$ and define the quotient slope surpluses by

```math
\epsilon_j=\kappa(F_j)-\kappa(E),
\qquad \kappa(F)=\frac{2\pi\mu(F)}{\mathrm{Vol}(X)}.
```

**Uniform relative asymptotics (Theorem 1.1).** Locally, the extension is
stable exactly when every quotient surplus is positive. Its first $k$
trace-free Dolbeault eigenvalues agree, to relative error
$O(\sum_{j=1}^k\epsilon_j)$, with the eigenvalues of

```math
S(\epsilon)=\mathrm{diag}(\epsilon_1,\ldots,\epsilon_k)
 +\frac{1}{r_0}vv^*,\qquad v_j=\sqrt{r_j\epsilon_j}.
```

All remaining eigenvalues stay uniformly positive. The relative estimate
is uniform across all positive gap ratios, including arbitrarily separated
rates of vanishing. The leading spectrum depends only on the ranks and
slopes. Exact divisibility of the scalar Hermitian–Einstein equations by
the individual extension weights makes this uniform estimate possible.

**The second-order term (Theorem 4.1).** For two factors, with positive
Einstein-constant difference $\delta$, the actual Hermitian–Einstein metrics
satisfy

```math
\lambda_1=\delta-D\delta^2+o(\delta^2),\qquad D\geq0.
```

The coefficient $D$ is given by the Green operator at the graded limit.
After rescaling the extension and making unitary identifications with the
graded limit, the normalized first eigensection is determined through order
$\delta$ in every spatial $C^m$ norm. This identifies the earlier harmonic-family
coefficient for equilibrium metrics with the slope gap as parameter.

**Toric examples (Section 5).** Tangent bundles of smooth toric Fano
varieties realize arbitrarily many independent small spectral scales.
The geometric extensions and local chambers come from
[Local stability chambers for toric tangent bundles](https://doi.org/10.5281/zenodo.23681809).
The joint limit in polarization and HYM flow time, repeated factors, and
extension graphs with cycles are outside the theorems.

## Build and finite computations

```bash
# From the repository root:
make paper3         # build/paper3/main.pdf; leaves the released PDF unchanged
make checks-paper3  # exact matrix computations
make integrity     # verify both distribution manifests
```

The five TeX files are self-contained, with an inline bibliography and no
external data. The [matrix checks](checks/README.md) use only the Python
standard library. They verify finite algebraic identities and examples;
the analytic estimates are proved in the manuscript.
The [source guide](paper/README.md) describes a direct LaTeX build.
The [SHA-256 manifest](SHA256.json) covers this distribution's sources,
released PDF, computation, guides, citation files and rights notice.

## Rights

Copyright 2026 Bernd Johannes Wuebben. All rights reserved for this paper,
its sources, supporting code and original prose. See [LICENSE](LICENSE).
The earlier flow paper retains its existing licenses.

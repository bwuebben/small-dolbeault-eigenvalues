# Small Dolbeault eigenvalues

Two research preprints by **Bernd Johannes Wuebben**
([ORCID](https://orcid.org/0009-0002-4588-3284)), with LaTeX sources and
supporting finite-dimensional computations.

Both papers study the Dolbeault Laplacian on trace-free endomorphisms of
holomorphic vector bundles. The flow paper follows a fixed bundle as time
tends to infinity. The stability-wall paper studies Hermitian–Einstein
equilibrium metrics as the Kähler class varies.

## Papers

| Manuscript | Current version | Read and cite |
|---|---|---|
| **[Small Dolbeault eigenvalues near stability walls](paper3/README.md)** | 27 October 2026 · 13 pages | [PDF](paper3/paper/main.pdf) · [Abstract](https://bwuebben.github.io/papers/small-dolbeault-stability-walls/) · [Zenodo:23717764](https://zenodo.org/records/23717764) |
| **[Small Dolbeault eigenvalues along Hermitian Yang–Mills flow](paper/README.md)** | 30 September 2026 · 45 pages; initial version 20 September | [PDF](paper/small-dolbeault-eigenvalues.pdf) · [Abstract](https://bwuebben.github.io/papers/small-dolbeault-eigenvalues/) · [Zenodo:23643223](https://zenodo.org/records/23643223) |

**Near stability walls.** For star extensions of pairwise nonisomorphic
factors that are stable of equal slope at the limiting Kähler class,
every small Hermitian–Einstein eigenvalue is determined to uniform relative
error by a matrix depending only on factor ranks and quotient slope
surpluses. The estimate holds across all positive gap ratios. The paper
also gives a two-factor second-order expansion and toric examples with
arbitrarily many independent spectral scales. Its
[guide](paper3/README.md) states the hypotheses and claim boundaries.
ArXiv submission is pending.

**Along Hermitian Yang–Mills flow.** A relative comparison with a weighted
graph Laplacian describes every small eigenvalue for simple semistable
curve bundles with distinct stable factors, and for two-step extensions
in higher dimensions. Under a positive-multiplier condition, the paper
determines every leading decay rate and coefficient, including which
coefficients retain initial-metric information. It also treats harmonic
families with repeated factors and proves sharp slope inequalities.
The [flow-paper guide](paper/README.md) gives the detailed results and attribution.

The papers concern different limiting problems. Neither establishes a
joint polarization/flow-time limit. The numbering `paper3/` is retained
for stable links; this repository contains the two public papers above.

## Repository map

| Location | Contents |
|---|---|
| [paper/](paper/) | Flow manuscript: released PDF, 13 TeX files and detailed guide |
| [checks/](checks/README.md) | Seven computations for the flow paper |
| [paper3/](paper3/README.md) | Stability-wall manuscript: released PDF, five TeX files, exact matrix computation, citations and rights |
| [CITATION.bib](CITATION.bib) | Separate BibTeX entries for both manuscripts |
| [CITATION.cff](CITATION.cff) | Repository metadata; preferred citation is the flow paper |
| [SHA256.json](SHA256.json), [paper3/SHA256.json](paper3/SHA256.json) | SHA-256 manifests for the two distributions |

The stability-wall examples use the geometric constructions in
*[Local stability chambers for toric tangent bundles](https://bwuebben.github.io/papers/toric-local-stability-chambers/)*
([source and supplement](https://github.com/bwuebben/toric-tangent-stability/tree/main/paper3),
[Zenodo:23682140](https://zenodo.org/records/23682140)).

## Build and verify

The computations require Python 3 and use only its standard library.
Building the PDFs requires `make` and a TeX Live installation with
`latexmk`, `amsart`, `geometry`, `hyperref` and `booktabs`. No external data is needed.
From the repository root:

```bash
make checks        # all eight computations across both papers
make integrity     # verify both SHA-256 manifests
make papers        # build both manuscripts
make verify        # computations, integrity and both builds
```

Individual targets are available:

```bash
make checks-flow   # seven flow computations
make checks-paper3 # exact matrix identities for stability walls
make paper         # flow PDF -> build/main.pdf
make paper3        # stability-wall PDF -> build/paper3/main.pdf
make clean         # remove build products
```

Builds leave the released PDFs unchanged. The computations verify finite
algebra and model examples; they do not verify the analytic proofs.
The [flow checks](checks/README.md) and
[stability-wall checks](paper3/checks/README.md) describe their coverage.
GitHub Actions runs the computations, checksum checks, citation validation
and both LaTeX builds on pushes and pull requests.

## Citation and rights

Cite the particular manuscript used, with its version DOI:

- Stability walls: [10.5281/zenodo.23717764](https://doi.org/10.5281/zenodo.23717764)
  · [BibTeX](paper3/CITATION.bib) · [CFF](paper3/CITATION.cff).
- HYM flow: [10.5281/zenodo.23643223](https://doi.org/10.5281/zenodo.23643223)
  · [BibTeX](CITATION.bib) · [CFF](CITATION.cff).

GitHub's **Cite this repository** uses the flow paper as its preferred
citation. Use the separate entry above when citing the stability-wall paper.
Manuscript version dates are distinct from archive deposit dates.

Rights differ between the papers. The flow manuscript and original prose
are CC BY 4.0; its supporting code is MIT licensed ([LICENSE](LICENSE)).
The stability-wall manuscript, sources, supporting code and original prose
are all rights reserved ([paper3/LICENSE](paper3/LICENSE)).

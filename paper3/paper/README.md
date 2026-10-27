# Manuscript

*Small Dolbeault eigenvalues near stability walls*, Bernd Johannes Wuebben,
27 October 2026, 13 pages. [Read the PDF](main.pdf).

From the repository root, run `make paper3` to build into
`build/paper3/main.pdf` without changing the released PDF. Alternatively,
from this directory:

```bash
mkdir -p ../../build/paper3
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../../build/paper3 main.tex
```

The main file includes `geometry.tex`, `spectrum.tex`, `two_factor.tex` and
`toric.tex`; its bibliography is inline. No external data or generated
inputs are required. The [paper guide](../README.md) states the results,
their scope, and the supporting finite computations.

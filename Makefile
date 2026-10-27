.DEFAULT_GOAL := help
.PHONY: help paper paper3 papers checks checks-flow checks-paper3 integrity verify clean

PYTHON ?= python3
FLOW_CHECKS := $(sort $(wildcard checks/check_*.py))

help:
	@printf '%s\n' \
	  'make paper          Build the flow manuscript into build/main.pdf' \
	  'make paper3         Build stability walls into build/paper3/main.pdf' \
	  'make papers         Build both manuscripts' \
	  'make checks         Run all eight computations across both papers' \
	  'make checks-flow    Run the seven flow computations' \
	  'make checks-paper3  Run the exact stability-wall matrix computation' \
	  'make integrity      Verify both SHA-256 manifests' \
	  'make verify         Run computations, integrity checks and both builds' \
	  'make clean          Remove build/'

paper:
	@command -v latexmk >/dev/null 2>&1 || \
	  { printf 'latexmk is required; install TeX Live before building.\n' >&2; exit 1; }
	@mkdir -p build
	cd paper && latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
	  -outdir=../build main.tex
	@printf 'PDF: build/main.pdf\n'

paper3:
	@command -v latexmk >/dev/null 2>&1 || \
	  { printf 'latexmk is required; install TeX Live before building.\n' >&2; exit 1; }
	@mkdir -p build/paper3
	cd paper3/paper && latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
	  -outdir=../../build/paper3 main.tex
	@printf 'PDF: build/paper3/main.pdf\n'

papers: paper paper3

checks: checks-flow checks-paper3

checks-flow:
	@set -e; for check in $(FLOW_CHECKS); do \
	  printf 'Running %s\n' "$$check"; $(PYTHON) "$$check"; done
	@printf '%s flow computations completed successfully.\n' "$(words $(FLOW_CHECKS))"

checks-paper3:
	$(PYTHON) paper3/checks/check_rank_matrix.py

integrity:
	$(PYTHON) -c 'import hashlib, json; from pathlib import Path; manifests = [Path("SHA256.json"), Path("paper3/SHA256.json")]; bad = [str(m.parent / p) for m in manifests for p, h in json.loads(m.read_text()).items() if hashlib.sha256((m.parent / p).read_bytes()).hexdigest() != h]; print("Checksum mismatches: " + ", ".join(bad) if bad else "Both SHA-256 manifests verified."); raise SystemExit(bool(bad))'

verify: checks integrity papers

clean:
	rm -rf build

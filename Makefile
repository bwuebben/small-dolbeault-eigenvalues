.DEFAULT_GOAL := help
.PHONY: help paper checks verify clean

PYTHON ?= python3
CHECKS := $(sort $(wildcard checks/check_*.py))

help:
	@printf '%s\n' \
	  'make paper    Build the manuscript into build/main.pdf' \
	  'make checks   Run the seven computational checks' \
	  'make verify   Run the checks, then build the manuscript' \
	  'make clean    Remove build/'

paper:
	@command -v latexmk >/dev/null 2>&1 || \
	  { printf 'latexmk is required; install TeX Live before building.\n' >&2; exit 1; }
	@mkdir -p build
	cd paper && latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
	  -outdir=../build main.tex
	@printf 'PDF: build/main.pdf\n'

checks:
	@set -e; for check in $(CHECKS); do \
	  printf 'Running %s\n' "$$check"; $(PYTHON) "$$check"; done
	@printf '%s checks completed successfully.\n' "$(words $(CHECKS))"

verify: checks paper

clean:
	rm -rf build

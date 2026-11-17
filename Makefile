PYTHON ?= python3

.PHONY: all verify pdf clean
all: verify pdf

verify:
	$(PYTHON) scripts/verify_results.py
	$(PYTHON) scripts/verify_strengthenings.py
	$(PYTHON) scripts/verify_partition_bounds.py

pdf:
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error optimal_coupon_packages.tex
	$(PYTHON) scripts/check_latex_log.py paper/optimal_coupon_packages.log

clean:
	cd paper && latexmk -c optimal_coupon_packages.tex

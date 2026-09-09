SPHINXOPTS ?=
SPHINXBUILD ?= sphinx-build
SOURCEDIR = .
BUILDDIR = _build

.PHONY: help strict linkcheck clean

help:
	@echo "HRIBO documentation targets:"
	@echo "  strict     rebuild HTML with warnings and missing references as errors"
	@echo "  linkcheck  rebuild the external-link inventory with warnings as errors"
	@echo "  clean      remove generated documentation output"

strict:
	@$(SPHINXBUILD) -E -W --keep-going -n -b html \
		"$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)

linkcheck:
	@$(SPHINXBUILD) -E -W --keep-going -n -b linkcheck \
		"$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS)

clean:
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

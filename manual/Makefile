LATEXMK ?= latexmk

SOURCE := manual.tex
ASSETS := ../HRIBO.png
JOBNAME := Manual_HRIBO
BUILD_DIR := build
PDF := $(JOBNAME).pdf

LATEX_FLAGS := -pdf \
	-interaction=nonstopmode \
	-halt-on-error \
	-file-line-error \
	-outdir=$(BUILD_DIR) \
	-jobname=$(JOBNAME)

.PHONY: all clean

all: $(PDF)

$(PDF): $(SOURCE) $(ASSETS) Makefile
	mkdir -p "$(BUILD_DIR)"
	$(LATEXMK) $(LATEX_FLAGS) "$(SOURCE)"
	cp "$(BUILD_DIR)/$(PDF)" "$(PDF)"

clean:
	$(LATEXMK) -C -outdir="$(BUILD_DIR)" -jobname="$(JOBNAME)" "$(SOURCE)"
	rm -rf -- "$(BUILD_DIR)"

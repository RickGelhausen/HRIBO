# HRIBO PDF manual

`manual.tex` is the editable source and `Manual_HRIBO.pdf` is its generated
user-facing PDF. Do not edit the PDF directly.

On Debian or Ubuntu, install a suitable local toolchain with:

```console
sudo apt-get install latexmk texlive-latex-extra texlive-fonts-recommended
```

Build from anywhere in the repository:

```console
make -C manual
```

Temporary LaTeX files are written to `manual/build/`. Remove them with
`make -C manual clean`.

The `Manual PDF` GitHub Actions workflow compiles the source for every relevant
pull request. On pushes to `development` or `master`, it also commits a changed
`Manual_HRIBO.pdf` back to the branch and publishes the compiled PDF as a
workflow artifact.

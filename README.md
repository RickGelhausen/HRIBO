<img src="HRIBO.png" width="620" alt="HRIBO">

# HRIBO

[![CI](https://github.com/RickGelhausen/HRIBO/actions/workflows/ci.yaml/badge.svg?branch=development)](https://github.com/RickGelhausen/HRIBO/actions/workflows/ci.yaml)
[![Documentation Status](https://readthedocs.org/projects/hribo/badge/?version=latest)](https://hribo.readthedocs.io/)
[![Snakemake](https://img.shields.io/badge/Snakemake-9.25.2-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE)

HRIBO (High-throughput annotation by Ribo-seq) is a reproducible Snakemake
workflow for bacterial ribosome-profiling data. It provides read processing and
quality control, strand-aware coverage tracks, metagene analysis and TIS advice,
feature counting, Reparation and optional DeepRibo ORF prediction, matched
RNA/Ribo differential analysis, and consolidated result tables.

> **Development status:** HRIBO 2.0 is under active validation. The automated
> suite and the production container boundaries are exercised in CI. A
> representative real-data comparison remains a release gate; see the
> [validation protocol](docs/real-data-validation.rst).

## Quick start

The tested path is Linux x86-64 with Conda or Micromamba and Apptainer. Keep the
workflow checkout separate from each analysis directory.

```console
git clone --branch development --single-branch \
  https://github.com/RickGelhausen/HRIBO.git /path/to/HRIBO
micromamba create --name hribo --file /path/to/HRIBO/environment.linux-64.pin.txt
micromamba activate hribo

mkdir -p /path/to/my-analysis/config
cp /path/to/HRIBO/config/config.yaml /path/to/my-analysis/config/
cp /path/to/HRIBO/config/samples.tsv /path/to/my-analysis/config/
cd /path/to/my-analysis
```

Edit the copied configuration and sample sheet, then add a matching genome
FASTA, GFF3/GTF annotation, and gzip-compressed FASTQ files. Check the complete
DAG before starting work:

```console
/path/to/HRIBO/run_hribo.sh \
  --configfile "$PWD/config/config.yaml" \
  --cores 20 \
  --dry-run all

/path/to/HRIBO/run_hribo.sh \
  --configfile "$PWD/config/config.yaml" \
  --cores 20 \
  --rerun-incomplete all
```

Select deliverables with `workflowSettings.stages`, or override them for one
run, for example `--config stages=mapping,tracks`. The `preprocessing` preset
selects trimming, mapping, and QC; `full` selects every stage and therefore
requires a valid matched differential-expression design.

For SLURM, activate the launcher environment, configure
`workflow/profiles/slurm/config.yaml` for the site, and invoke
`/path/to/HRIBO/slurm_run.sh` from the analysis directory.

## Documentation

The maintained documentation source and build configuration now live in
[`docs/`](docs/index.rst); the former `HRIBO_ReadTheDocs` repository is no
longer a documentation source. Its complete history is preserved on the
`archive/hribo-readthedocs` branch. The hosted site may continue to show the
legacy version until a project administrator completes the
[Read the Docs cutover](docs/development.rst#read-the-docs-cutover). Start with:

- [installation and execution](docs/getting-started.rst)
- [sample-sheet format](docs/samples.rst)
- [configuration](docs/configuration.rst)
- [stage selection](docs/stages.rst)
- [output catalogue](docs/outputs.rst)
- [result-table reference](docs/table-reference.rst)
- [historical public example data](docs/historical-example-data.rst)
- [migration from HRIBO 1.8](docs/migration-1.8-to-2.0.rst)

## Citation

If HRIBO contributes to published work, please cite:

> Gelhausen R, Svensson SL, Froschauer K, et al. HRIBO: high-throughput
> analysis of bacterial ribosome profiling data. *Bioinformatics*.
> 2021;37(14):2061–2063.
> <https://doi.org/10.1093/bioinformatics/btaa959>

Machine-readable metadata is provided in [`CITATION.cff`](CITATION.cff).
HRIBO is licensed under GPL-3.0; see [`LICENSE`](LICENSE).

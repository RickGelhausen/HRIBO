"""
Which parts of the pipeline a run should produce.

The config names either a preset ("full", "mapping", ...) or an explicit list of
stages. A stage only says what is *requested*; Snakemake still builds everything
a request depends on, so asking for "predictions" alone trims, filters, maps and
counts reads on the way there. Conversely, deselecting a stage removes its
outputs from the DAG, which is how a run can stop at the BAM files.

Only the standard library is used, so that this module can be imported from the
Snakefile at parse time and from the preflight checks.

Author: Rick Gelhausen
"""

from __future__ import annotations

# Stage name -> what asking for it puts in the output directory. The order is
# roughly the order of the pipeline and is preserved when the stages are
# resolved, so error messages and the startup banner read sensibly.
STAGE_DESCRIPTIONS = {
    "trimming": "adapter removal, plus FastQC of the raw and trimmed reads",
    "mapping": "the final alignments in maplink/ with their indices",
    "qc": "the MultiQC report over every quality control step",
    "tracks": "bigwig coverage tracks",
    "genome_tracks": "start, stop, alternative start codon and RBS GFF tracks",
    "readcounts": "read count and annotation spreadsheets",
    "metagene": "metagene profiles and read length statistics",
    "tis_advisor": "read length and P-site offset advice for a TIS caller",
    "correlation": "the Spearman correlation heatmap between libraries",
    "pca": "the PCA plot of the read counts",
    "predictions": "annotation independent ORF predictions",
    "differential_expression": "xtail, riborex and deltaTE contrast tables",
    "overview": "the combined overview spreadsheet",
}

STAGE_NAMES = list(STAGE_DESCRIPTIONS)

# Input requirements include the dependencies Snakemake builds on the way to a
# requested stage.  Trimming consumes reads but no reference; genome motif
# tracks consume the genome but neither reads nor an annotation.  Every other
# stage descends from mapping and therefore needs all three input classes.
FASTQ_REQUIRED_STAGES = frozenset(STAGE_NAMES) - {"genome_tracks"}
GENOME_REQUIRED_STAGES = frozenset(STAGE_NAMES) - {"trimming"}
ANNOTATION_REQUIRED_STAGES = frozenset(STAGE_NAMES) - {
    "trimming",
    "genome_tracks",
}

# Stages that need Ribo-seq libraries. Without a RIBO library there is nothing
# for the ORF predictors to work on, and everything built on top of them.
RIBO_ONLY_STAGES = {"predictions", "differential_expression", "overview"}

# Metagene profiles can also be built from translation-initiation or
# translation-termination profiling.  Keep this separate from RIBO_ONLY_STAGES:
# a TIS-only project has useful metagene evidence even though the general RIBO
# prediction and differential-expression stages are unavailable.
RIBO_LIKE_METHODS = {"RIBO", "TIS", "TTS"}
RIBO_LIKE_STAGES = {"metagene", "tis_advisor"}

# Named combinations, so that the common cases stay a single word. Preset names
# deliberately do not collide with stage names; a stage is always looked up
# first, so "mapping" means the stage and nothing more than the stage.
PRESETS = {
    "full": STAGE_NAMES,
    "preprocessing": ["trimming", "mapping", "qc"],
}


class StageError(ValueError):
    """Raised when the requested stages cannot be interpreted."""


def _as_list(value) -> list[str]:
    """Accept a preset name, a comma-separated string or a list of stages.

    The comma-separated form is what reaches us from the command line, where
    `--config stages=mapping,tracks` overrides the config file for a single run.
    """
    if isinstance(value, str):
        return [part.strip() for part in value.split(",") if part.strip()]
    if isinstance(value, (list, tuple)):
        return [str(part).strip() for part in value if str(part).strip()]
    raise StageError(
        f"workflowSettings.stages must be a preset name or a list of stages, got {value!r}"
    )


def resolve_stages(
    config,
    has_ribo: bool = True,
    has_ribo_like: bool | None = None,
) -> list[str]:
    """The stages to build, in pipeline order and with the presets expanded.

    A top level `stages` key wins over `workflowSettings.stages`, which is how
    `snakemake --config stages=mapping` overrides the config file without
    editing it.

    `has_ribo` drops stages that specifically require a RIBO library.
    `has_ribo_like` controls metagene stages, which can use RIBO, TIS or TTS.
    It defaults to `has_ribo` for backwards-compatible callers that only know
    whether a project has ordinary Ribo-seq.
    """
    requested = config.get("stages")
    if requested is None:
        requested = config.get("workflowSettings", {}).get("stages", "full")

    names = _as_list(requested)
    if not names:
        raise StageError("No stages requested. Set workflowSettings.stages to a preset or a list of stages.")

    selected = set()
    for name in names:
        if name in STAGE_DESCRIPTIONS:
            selected.add(name)
        elif name in PRESETS:
            selected.update(PRESETS[name])
        else:
            known = ", ".join(sorted(set(PRESETS) | set(STAGE_DESCRIPTIONS)))
            raise StageError(f"Unknown stage or preset {name!r}. Known values: {known}")

    if has_ribo_like is None:
        has_ribo_like = has_ribo

    if not has_ribo:
        selected -= RIBO_ONLY_STAGES
    if not has_ribo_like:
        selected -= RIBO_LIKE_STAGES

    return [name for name in STAGE_NAMES if name in selected]

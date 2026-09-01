# HRIBO remediation progress

**Started:** 2026-09-01

**Source audit:** [PROJECT_AUDIT.md](PROJECT_AUDIT.md)

**Target:** a trustworthy, reproducible HRIBO 2.0 release

This file is the working record for changes made in response to the project audit. It is intentionally separate from the audit: the audit records the state that was found, while this document records remediation decisions, code changes, tests, and remaining work.

The repository already contained an unfinished stage-selection overhaul when remediation began. Those pre-existing changes are being preserved and are not attributed to this remediation unless a row below explicitly says otherwise.

## Status legend

| Status | Meaning |
|---|---|
| Open | Confirmed issue; no implementation started. |
| In progress | Being implemented or reviewed. |
| Fixed | Code and focused regression coverage are complete. |
| Verified | Focused tests and the full available verification suite pass. |
| Deferred | Intentionally postponed, with a reason recorded. |
| Removed | Unsupported/dead functionality was removed instead of repaired. |

## Current remediation batch

No correctness implementation is active at this checkpoint. Batches 1–4 are verified. The remaining open work is release engineering, dependency reproducibility, CI/integration execution, legacy-path decisions, and documentation; those items remain explicit in the ledger rather than being presented as complete.

## Completed remediation batches

| Batch | IDs | Result |
|---|---|---|
| 1 — scientific correctness and failure signaling | META-001, DIFF-001, OVERVIEW-001 | Verified on 2026-09-01. |
| 2 — input, annotation, and metagene edge cases | PATH-001, REF-001, ANNOT-001, STAGE-001, META-002, META-004 | Verified on 2026-09-01. |
| 3 — workflow invocation, configuration, validation, and small designs | RUN-001, QC-001, SCRIPT-001, CONFIG-001, ANNOT-002, CONTRAST-001, DIFF-002, DIFF-003, PCA-001, META-003, VALID-001 | Verified on 2026-09-01. |
| 4 — prediction validity, determinism, and zero-result propagation | PRED-001, PRED-002, REPAIR-001 | Verified on 2026-09-01. |

## Remediation ledger

### Scientific correctness and failure signaling

| ID | Priority | Finding | Status | Notes |
|---|---:|---|---|---|
| DIFF-001 | Critical | deltaTE can report success after the underlying analysis fails. | Verified | Uses a checked runner; tool exits and missing outputs fail, stale artifacts are cleared, and under-replicated designs fail preflight. |
| OVERVIEW-001 | High | Explicit differential-expression contrasts are ignored in the overview. | Verified | Explicit order/orientation is preserved and populated values are checked. |
| META-001 | High | Stop profiles are misoriented; global minus-start and plus-stop coverage can be dropped. | Verified | One strand transform is shared by both anchors; 16 asymmetric mapping cases pass. |
| ANNOT-001 | High | GTF conversion has undefined/stale variables, unsafe format handling, missing-output success, and incorrect child coordinates. | Verified | Replaced with a deterministic, atomic converter that preserves feature bounds, reports unsupported/mixed input, and emits strictly valid GFF3. |
| DIFF-002 | Medium | Differential-expression inputs do not validate matched RIBO/RNA sample widths. | Verified | Validation is contrast-scoped and rejects unequal selected RIBO/RNA replicate counts before R tools run. |
| DIFF-003 | Medium | Self-contrasts create duplicated count matrices and a meaningless one-level statistical design. | Verified | Preflight rejects contrasts whose two condition names are identical. |
| PCA-001 | Medium | PCA unconditionally requests three components. | Verified | Uses the effective numerical rank, supports one- and two-component plots, and requires at least two libraries only when PCA is selected. |
| META-002 | Medium | Empty start/stop evidence can call `min()` on an empty list; retained annotation keys can also lack a read-interval entry. | Verified | Missing read indexes become zero evidence; empty/start-only/stop-only inputs remain explicit and window-normalized zero columns stay finite. |
| META-003 | Medium | `lengthCutoff`, PDF output, and multi-color handling do not behave as configured. | Verified | The cutoff reaches filtering, PDF is supported, and configured colors control per-length traces with deterministic cycling/fallback. |
| META-004 | Medium | Metagene boundary filtering retains strand/edge windows that extend outside a contig. | Verified | Explicit inclusive start/stop windows now use strand-aware asymmetric geometry and exact zero-based contig bounds. |

### Workflow reliability and edge cases

| ID | Priority | Finding | Status | Notes |
|---|---:|---|---|---|
| RUN-001 | High | `tis_advisor.py` is invoked directly despite not being executable in Git. | Verified | The rule invokes the configured Python interpreter explicitly, quotes arguments, and captures the real command log. |
| PATH-001 | High | Absolute FASTQ paths are prefixed with the working directory and become broken link targets. | Verified | A replacement-safe stager resolves the source once and creates correct absolute links for relative, absolute, and space-containing paths. |
| REF-001 | High | Gzipped references pass validation but are copied without decompression. | Verified | Reference staging detects gzip by magic bytes and atomically materializes plain downstream FASTA/GFF files. |
| STAGE-001 | High | RNA-only metagene selection schedules an empty `-a` argument list. | Verified | The resolver removes metagene/TIS-advisor targets when no RIBO-like assay exists while retaining valid TIS/TTS-only runs. |
| PRED-001 | High | Empty DeepRibo results fail to create every declared output. | Verified | Both GFF outputs are unconditional, and schema-bearing empty results propagate through counting, mapping, and Excel output. |
| PRED-002 | High | Append-mode and undeclared intermediate files can preserve stale predictions. | Verified | Aggregators replace outputs atomically/deterministically, DeepRibo updates its pair transactionally, and every prediction CDS is strict GFF3 with phase `0`. |
| QC-001 | Medium | Paired trimmed FastQC passes both mates to both commands. | Verified | Each named mate is processed exactly once, with quoted paths, declared outputs, and captured logs. |
| REPAIR-001 | Medium | Reparation aggregation order depends on a Python set and lacks a final sort. | Verified | A deterministic maximum-probability winner supplies coherent metadata; evidence and final coordinates are sorted. |
| SCRIPT-001 | Medium | `enrich_annotation.py` calls `sys.exit()` without importing `sys`. | Verified | The error path imports `sys`, exits nonzero, and reports the malformed identifier. |
| LEGACY-001 | Medium | `annotationBigBed`, `colorBigWig`, and several logging/provenance paths are broken or orphaned. | In progress | Metagene/TIS/FastQC logging was repaired; PEAR/Reparation capture, contrast-marker provenance, and the two orphan visualization rules remain. |
| STAGE-002 | Medium | The trimming stage targets QC pages but does not retain trimmed/assembled reads as stage deliverables. | Open | Define whether this is a trimming-output or trimming-QC stage. |
| OUTPUT-001 | Medium | Combined updated annotation rules exist but no stage exposes their output. | Open | Expose the updated annotation under predictions or remove the dormant path. |
| DEPEND-001 | Medium | Most directly invoked scripts are not declared as Snakemake job dependencies. | Open | Prefer `script:` or explicit `workflow.source_path(...)` inputs. |
| PORT-001 | Medium | Final `maplink/` BAM symlinks are absolute and break when a result directory moves. | Open | Replace them with portable relative links or real final files. |

### Validation and configuration integrity

| ID | Priority | Finding | Status | Notes |
|---|---:|---|---|---|
| CONFIG-001 | Medium | Configured alternative start codons are ignored. | Verified | The motif rule uses the configured codon list; matching is literal, case-insensitive, overlapping, and safe for an empty list. |
| ANNOT-002 | Medium | Malformed annotation rows are collected but never reported. | Verified | Parsing returns malformed-row details and preflight emits a blocking `ANNOTATION_MALFORMED_ROWS` diagnostic. |
| CONTRAST-001 | Medium | Auto-contrasts include conditions from assays unsupported by the differential tools. | Verified | Automatic contrasts use the deterministic intersection of RIBO and RNA condition sets. |
| VALID-001 | Medium | Validation is not stage-aware. | Verified | File preflight now follows actual stage dependencies: trimming skips references, genome tracks skip annotation/FASTQ, and mapping/full still require all inputs. |

### Reproducibility, packaging, and maintenance

| ID | Priority | Finding | Status | Notes |
|---|---:|---|---|---|
| REPRO-001 | High | Containers, the DeepRibo model, and Swiss-Prot inputs are mutable. | Open | Pin digests/releases/commits and record checksums. |
| ENV-001 | High | Old and incompletely declared environments weaken clean installation. | Open | Solve and pin every supported environment after code-path decisions. |
| SHELL-001 | Medium | Several user-controlled paths are interpolated without shell quoting. | In progress | New and repaired paths are quoted or use Python staging; audit the remaining legacy rules repository-wide. |
| TEST-001 | High | CI lacks a real end-to-end workflow fixture. | Open | Add minimal partial and full execution fixtures. |
| TEST-002 | Medium | Optional test dependencies can silently skip major spreadsheet/TIS coverage. | Open | Provide a complete test environment and assert the expected collected/skipped counts in CI. |
| LINT-001 | Medium | Repository-wide Ruff debt remains outside the clean changed-file baseline. | Open | Fix or explicitly baseline all outstanding findings, then enforce it in CI. |
| PACKAGE-001 | Medium | Wheel/package discovery points at nonexistent directories. | Open | Decide workflow-only versus installable package first. |
| DOCS-001 | Medium | README, manual, ReadTheDocs, tags, and version metadata disagree. | Open | Rebuild documentation after supported behavior stabilizes. |

## Change log

| Date | IDs | Files | Change | Evidence |
|---|---|---|---|---|
| 2026-09-01 | — | `REMEDIATION_PROGRESS.md` | Created the remediation ledger and first prioritized batch. | Ledger reviewed against `PROJECT_AUDIT.md`. |
| 2026-09-01 | META-001 | `workflow/scripts/lib/metagene.py`, `workflow/scripts/tis_advisor.py`, `tests/test_metagene_lib.py` | Unified plus-direct/minus-reverse transcript coordinates, ordered global slices, and gave stop profiles their correct axis. | 16 anchor/strand/mode cases plus TIS end-to-end tests; independent semantic review approved. |
| 2026-09-01 | OVERVIEW-001 | `workflow/scripts/generate_excel_overview.py`, `tests/test_excel_outputs.py`, `tests/golden/overview.csv` | Preserved explicit contrast order/orientation and retained deterministic inference only when `-c` is omitted. | `B-A` headers and populated values tested; non-contrast golden cells unchanged. |
| 2026-09-01 | DIFF-001 | `workflow/rules/diffex_deltate.smk`, `workflow/scripts/run_deltate.sh`, `workflow/scripts/lib/checks.py`, `tests/test_deltate_rule.py`, `tests/test_validation.py` | Removed failure suppression and fake outputs, added checked execution, stale cleanup, replication preflight, and contrast-scoped validation. | Fake `DTEG.R` exercises exit 23, missing output, stale output, header-only success, and paths with spaces; independent review approved. |
| 2026-09-01 | PATH-001, REF-001 | `workflow/scripts/stage_input.py`, `workflow/rules/trimming.smk`, `workflow/rules/preprocessing.smk`, `tests/test_input_staging.py` | Centralized replacement-safe FASTQ links and atomic plain-text reference materialization, including gzip magic-byte detection. | Relative, absolute, spaced, rerun, dangling-link, plain, gzip, and corrupt-gzip cases covered. |
| 2026-09-01 | ANNOT-001 | `workflow/scripts/gtf2gff3.py`, `tests/test_gtf2gff3.py`, `tests/golden_gff/gtf2gff3.gff` | Rebuilt GTF conversion around explicit parentage, exact child bounds, valid phases/escaping, deterministic order, actionable conflicts, and atomic replacement. | 24 focused cases pass; golden output passes `gt gff3validator`. |
| 2026-09-01 | STAGE-001, META-002 | `workflow/scripts/lib/stages.py`, `workflow/scripts/lib/metagene.py`, `workflow/scripts/lib/misc.py`, `workflow/scripts/metagene_profiling.py`, `tests/test_stages.py`, `tests/test_metagene_lib.py`, `tests/test_metagene_settings.py` | Removed biologically inapplicable metagene targets and made missing/one-sided/empty evidence a supported result. | RNA-only, TIS/TTS-only, absent contig/strand, start-only, stop-only, and fully empty fixtures pass. |
| 2026-09-01 | RUN-001, META-003, QC-001 | `workflow/rules/metageneprofiling.smk`, `workflow/rules/qc.smk`, `workflow/scripts/lib/annotation.py`, `workflow/scripts/lib/io.py`, `workflow/scripts/metagene_profiling.py`, `workflow/scripts/tis_advisor.py`, `tests/test_metagene_settings.py`, `tests/test_qc_rules.py` | Corrected interpreter invocation, command logging, filtering arguments, cutoff propagation, output-format handling, multi-color argv construction, and paired FastQC inputs. | Rendered-rule assertions and functional setting tests pass; independent review found no remaining blocker. |
| 2026-09-01 | CONFIG-001, SCRIPT-001, ANNOT-002 | `workflow/rules/visualization.smk`, `workflow/scripts/motif_to_gff.py`, `workflow/scripts/enrich_annotation.py`, `workflow/scripts/lib/validation.py`, `workflow/scripts/lib/checks.py`, `tests/test_script_fixes.py`, `tests/test_validation.py` | Honored configured alternative starts, repaired enrichment diagnostics, and made malformed annotation rows block preflight. | Literal/overlapping/empty motif cases and malformed-annotation diagnostics pass. |
| 2026-09-01 | CONTRAST-001, DIFF-002 | `workflow/Snakefile`, `workflow/scripts/lib/checks.py`, `tests/test_validation.py`, `tests/test_stages.py` | Limited automatic contrasts to matched assay domains and rejected unequal selected replicate widths. | Mixed-assay and unequal-replicate fixtures pass. |
| 2026-09-01 | DIFF-003 | `workflow/scripts/lib/checks.py`, `tests/test_validation.py` | Rejected explicit self-contrasts before they can duplicate identical samples into a one-level design. | Focused validation suite passes; independent review reproduced the old failure and approved the guard. |
| 2026-09-01 | PCA-001 | `workflow/scripts/analyse_variance.R`, `workflow/scripts/plot_PCA.py`, `workflow/scripts/lib/checks.py`, `tests/test_pca.py`, `tests/test_validation.py` | Selected PCA components by effective numerical rank, added truthful 1D/2D fallback plots and zero-range padding, parsed correlation labels correctly, and added stage-scoped minimum-library validation. | 43 focused PCA/validation tests, R parsing, and targeted Ruff pass. |
| 2026-09-01 | VALID-001 | `workflow/scripts/lib/stages.py`, `workflow/scripts/validate.py`, `tests/test_validation.py` | Declared input classes required by each requested stage and gated filesystem preflight without weakening global configuration consistency checks. | Trimming/genome-track omission boundaries and required-file counter-cases pass within a 67-test validation/stage suite; independent dependency review approved. |
| 2026-09-01 | META-002, META-003, META-004 | `workflow/scripts/lib/annotation.py`, `workflow/scripts/lib/misc.py`, `workflow/scripts/lib/plotting.py`, `workflow/scripts/metagene_profiling.py`, `tests/test_metagene_lib.py`, `tests/test_metagene_settings.py` | Kept one-sided window normalization finite, corrected all strand/contig edge windows, and made custom colors affect actual Plotly traces. | 118 metagene/P-site/TIS tests pass, including raw/window, plus/minus, exact-boundary, asymmetric-window, and color-cycle cases. |
| 2026-09-01 | PRED-001, PRED-002, REPAIR-001 | `workflow/rules/deepribo.smk`, `workflow/rules/merge.smk`, prediction GFF/merge scripts, `call_featurecounts.py`, `map_reads_to_annotation.py`, `excel_utils.py`, prediction fixtures/goldens/tests | Removed append-state, made the DeepRibo output pair transactional, propagated schema-bearing empty results through downstream workbooks, normalized prediction attributes, moved DeepRibo rank/distance metadata out of phase, and made Reparation select one coherent winner. | 89 focused tests pass; all aggregate goldens contain CDS phase `0` and pass GenomeTools; rollback and independent semantic reviews pass. |

The prediction phase correction follows the [Sequence Ontology GFF3 1.26 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), which requires a `0`, `1`, or `2` phase for every CDS feature.

## Verification log

The audit baseline before remediation was:

- `pytest -q`: 213 passed.
- Python byte compilation: passed.
- YAML parsing: passed.
- `bash -n` for launch scripts: passed.
- Snakemake 9.25.2 partial, default, and full DAG dry-runs: passed.
- Ruff: 75 findings.
- Wheel build: failed because package discovery references nonexistent directories.

New verification runs will be recorded here with the exact scope and result. A finding is not promoted to **Verified** solely because a focused unit test passes; the full available suite must also remain green.

### 2026-09-01 — batch 1

- Focused affected suites: 93 passed.
- Full suite: 238 passed in 72.39 s.
- Snakemake 9.25.2 mapping dry-run: passed.
- Snakemake 9.25.2 full dry-run with eight matched RIBO/RNA libraries: passed.
- Python byte compilation: passed.
- YAML parsing: 29 files passed.
- Bash syntax: launch scripts and `run_deltate.sh` passed.
- Ruff on every changed Python file: passed.
- `git diff --check`: passed.

### 2026-09-01 — batches 2 and 3, focused checkpoint

- Combined non-prediction regression suites: 187 passed in 74.78 s.
- PCA and validation boundary suite: 43 passed in 6.57 s.
- GTF converter suite: 24 passed; strict GFF3 validation passed.
- Input-staging suite: 10 passed, including failure-safe reruns.
- Targeted Ruff, R syntax parsing, Python byte compilation, and `git diff --check`: passed at their respective implementation checkpoints.
- At this intermediate checkpoint, full-suite and Snakemake verification were pending; both are recorded as passed in the integrated checkpoint below.

### 2026-09-01 — integrated checkpoint, batches 1–4

- Full Python suite: **349 passed in 139.59 s**.
- Prediction/GFF/Excel suite: **89 passed in 77.66 s**.
- Metagene/P-site/TIS suite: **118 passed in 35.23 s**; an expanded independent review suite passed 132 tests.
- Validation/stage suite: **67 passed**; independent DAG-dependency review found no blocker.
- Snakemake **9.25.2** mapping dry-run: passed against the disposable eight-library fixture.
- Snakemake **9.25.2** full dry-run: passed for all stages, including differential expression and predictions.
- Ruff on every modified or untracked Python file: passed.
- Python byte compilation for `workflow/scripts` and `tests`: passed.
- YAML parsing: **29 files passed**.
- R syntax parsing for `analyse_variance.R`: passed.
- Bash syntax for both launchers and `run_deltate.sh`: passed.
- GenomeTools accepted the DeepRibo main/plus, Reparation, and converted-GTF goldens as valid GFF3; it emitted only nonblocking missing-`##sequence-region` warnings.
- `git diff --check`: passed.
- Independent reviews of metagene semantics, stage-aware validation, differential contrasts, and prediction aggregation: no blockers.

## Working principles

1. Correct scientific results and honest failure reporting come before refactoring or presentation work.
2. Every fixed defect receives regression coverage that would fail against the pre-fix implementation.
3. Empty/zero-result datasets are valid outcomes unless a tool explicitly requires otherwise.
4. Existing user overhaul changes are preserved and reviewed for compatibility rather than overwritten.
5. External tools, databases, models, and containers must eventually be versioned and checksum-verifiable.
6. Unsupported legacy paths will be removed explicitly rather than left present but broken.

## Next planned batch

Following this verified correctness checkpoint:

1. Decide which orphan legacy visualization and updated-annotation paths are supported versus removed (`LEGACY-001`, `OUTPUT-001`).
2. Harden remaining shell/dependency boundaries and portable final-output links (`SHELL-001`, `DEPEND-001`, `PORT-001`).
3. Add fixture-backed workflow execution and continuous integration (`TEST-001`, `TEST-002`, `LINT-001`).
4. Pin external artifacts/environments, settle packaging, and rebuild release documentation (`REPRO-001`, `ENV-001`, `PACKAGE-001`, `DOCS-001`).

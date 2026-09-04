"""Regression tests for shell quoting and portable final alignment links."""

import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import call_featurecounts
from lib.cli import command_option_values


REPO = Path(__file__).resolve().parent.parent
RULES = REPO / "workflow" / "rules"
STAGER = REPO / "workflow" / "scripts" / "stage_input.py"


def run_stager(operation, source, destination):
    return subprocess.run(
        [sys.executable, str(STAGER), operation, str(source), str(destination)],
        capture_output=True,
        text=True,
    )


def rule_body(path, name):
    text = path.read_text()
    start = text.index(f"rule {name}:")
    next_rule = text.find("\nrule ", start + 1)
    return text[start:] if next_rule == -1 else text[start:next_rule]


def shell_blocks(text):
    """Yield indented bodies belonging to Snakemake ``shell:`` directives."""

    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "shell:":
            continue
        indentation = len(line) - len(line.lstrip())
        body = []
        for candidate in lines[index + 1 :]:
            if candidate.strip():
                candidate_indent = len(candidate) - len(candidate.lstrip())
                if candidate_indent <= indentation:
                    break
            body.append(candidate)
        yield "\n".join(body)


def test_portable_link_survives_spaces_and_moving_the_result_tree(tmp_path):
    original = tmp_path / "original result with spaces"
    source = original / "bam" / "RIBO-A-1.bam"
    destination = original / "maplink" / "RIBO-A-1.bam"
    index = original / "maplink" / "RIBO-A-1.bam.bai"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"bam fixture")

    result = run_stager("portable-link", source, destination)
    assert result.returncode == 0, result.stderr
    index.write_bytes(b"index fixture")

    target = os.readlink(destination)
    assert not os.path.isabs(target)
    assert target == os.path.relpath(source, destination.parent)
    assert destination.read_bytes() == b"bam fixture"

    moved = tmp_path / "moved result with spaces"
    original.rename(moved)
    assert (moved / "maplink" / destination.name).read_bytes() == b"bam fixture"
    assert (moved / "maplink" / index.name).read_bytes() == b"index fixture"


def test_portable_link_rerun_replaces_an_existing_or_dangling_link(tmp_path):
    result_tree = tmp_path / "result"
    first = result_tree / "bam" / "first.bam"
    second = result_tree / "bam" / "second.bam"
    destination = result_tree / "maplink" / "sample.bam"
    first.parent.mkdir(parents=True)
    first.write_bytes(b"first")
    second.write_bytes(b"second")

    assert run_stager("portable-link", first, destination).returncode == 0
    first.unlink()
    result = run_stager("portable-link", second, destination)

    assert result.returncode == 0, result.stderr
    assert destination.is_symlink()
    assert destination.read_bytes() == b"second"
    assert os.readlink(destination) == os.path.relpath(second, destination.parent)


def test_maplink_rule_delegates_to_the_quoted_portable_link_helper():
    mapping = RULES / "mapping.smk"
    rule = rule_body(mapping, "maplink")

    assert "os.getcwd()" not in rule
    assert "ln -s" not in rule
    assert 'stager=str(SCRIPTS / "stage_input.py")' in rule
    assert (
        "python3 {input.stager:q} portable-link {input.bam:q} {output.bam:q}"
        in rule
    )

    bamindex = rule_body(RULES / "visualization.smk", "bamindex")
    assert '"maplink/{method}-{condition}-{replicate}.bam.bai"' in bamindex
    assert "samtools index -@ {threads} {input.bam:q}" in bamindex


def test_pear_and_reparation_capture_complete_logs_with_quoted_paths():
    pear = rule_body(RULES / "trimming.smk", "merge_fastq")
    assert pear.index("exec > {log:q} 2>&1") < pear.index("pear -n 10")
    for placeholder in (
        "{input.fastq1:q}",
        "{input.fastq2:q}",
        "{params.prefix:q}",
        "{params.assembled:q}",
        "{output.fastq:q}",
    ):
        assert placeholder in pear

    reparation = rule_body(RULES / "reparation.smk", "reparation")
    assert reparation.index("exec > {log:q} 2>&1") < reparation.index(
        "reparation.pl"
    )
    for placeholder in (
        "{input.bam:q}",
        "{input.genome:q}",
        "{input.gtf:q}",
        "{input.db:q}",
        "{params.prefix:q}",
    ):
        assert placeholder in reparation


def test_each_contrast_marker_writes_only_its_declared_output():
    rule = rule_body(RULES / "diffex_contrast.smk", "contrastInput")

    assert "for f in CONTRASTS" not in rule
    assert "run:" not in rule
    assert '"touch {output:q}"' in rule


def test_comma_separated_adapter_options_are_trimmed_and_empty_safe():
    assert command_option_values("-a", " AAA,CCC, , TTT ") == [
        "-a",
        "AAA",
        "-a",
        "CCC",
        "-a",
        "TTT",
    ]
    assert command_option_values("-g", "   ") == []


def test_featurecounts_keeps_spaced_paths_and_features_as_single_argv_tokens(
    tmp_path, monkeypatch
):
    inputs = tmp_path / "inputs with spaces"
    inputs.mkdir()
    annotation = inputs / "annotation with spaces.gff"
    bam = inputs / "reads with spaces.bam"
    output = tmp_path / "counts with spaces.raw"
    annotation.write_text(
        "chr1\ttest\tcustom feature\t1\t3\t.\t+\t0\tID=feature1\n"
    )
    bam.write_bytes(b"fake bam")
    calls = []

    def fake_featurecounts(command):
        calls.append(command)
        temporary = Path(command[command.index("-o") + 1])
        temporary.write_text(
            "# Program:featureCounts\n"
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\treads\n"
            "feature1\tchr1\t1\t3\t+\t3\t7\n"
        )
        return 0

    monkeypatch.setattr(call_featurecounts.subprocess, "call", fake_featurecounts)
    args = SimpleNamespace(
        bamfiles=[str(bam)],
        annotation=str(annotation),
        output=str(output),
        features=["custom feature"],
        strandness=1,
        threads=2,
        assign_to_all=True,
        assign_multi_mappers=True,
        with_fraction=True,
        diff_expr=False,
    )

    call_featurecounts.call_featureCounts(args)

    assert len(calls) == 1
    command = calls[0]
    assert command[command.index("-t") + 1] == "custom feature"
    assert command[command.index("-a") + 1] == str(annotation)
    assert command[-1] == str(bam)
    assert "with" not in command
    assert output.is_file()


def test_unsupported_browser_export_rules_and_environments_are_removed():
    visualization = (RULES / "visualization.smk").read_text()
    for name in (
        "annotationBed",
        "annotationBed6",
        "annotationBigBed",
        "colorBigWig",
        "colorGFF",
    ):
        assert f"rule {name}:" not in visualization

    for legacy_output in (
        "tracks/annotation.bb",
        "tracks/annotation-woGenes.gtf",
        "tracks/color/",
    ):
        assert legacy_output not in visualization

    environments = REPO / "workflow" / "envs"
    assert not (environments / "bed.yaml").exists()
    assert not (environments / "color.yaml").exists()


def test_genome_index_rule_owns_and_writes_its_declared_log():
    visualization = RULES / "visualization.smk"
    index_rule = rule_body(visualization, "genomeSamToolsIndex")
    size_rule = rule_body(visualization, "genomeSize")

    assert '"logs/genomeSamToolsIndex.log"' in index_rule
    assert "samtools faidx {input.genome:q} 2> {log:q}" in index_rule
    assert "genomeSamToolsIndex.log" not in size_rule
    assert "{log" not in size_rule


def test_plot_correlation_declares_matrix_path_and_captures_its_log():
    correlation = rule_body(RULES / "visualization.smk", "plotCorrelation")

    assert 'matrix="figures/SpearmanCorr_readCounts.tab"' in correlation
    assert "--outFileCorMatrix {output.matrix:q}" in correlation
    assert "--outFileCorMatrix SpearmanCorr_readCounts.tab" not in correlation
    assert '"logs/plotCorrelation.log"' in correlation
    assert "> {log:q} 2>&1" in correlation


def test_samuniq_rule_retains_strict_failure_propagation():
    mapping = rule_body(RULES / "mapping.smk", "samuniq")
    assert "set +e" not in mapping
    assert ".mapped" not in mapping
    assert "grep" not in mapping
    assert 'NH:i:1' in mapping


def test_samuniq_filter_keeps_only_mapped_unique_alignments_and_rejects_none(
    tmp_path,
):
    mapping = rule_body(RULES / "mapping.smk", "samuniq")
    command_line = next(
        line.strip()
        for line in mapping.splitlines()
        if line.strip().startswith("awk ") and "NH:i:1" in line
    )
    program = command_line.split("'", 2)[1].replace("{{", "{").replace("}}", "}")
    alignments = tmp_path / "alignments.sam"
    alignments.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "forward\t0\tchr1\t1\t255\t3M\t*\t0\t0\tAAA\tIII\tNH:i:1\n"
        "reverse\t16\tchr1\t2\t255\t3M\t*\t0\t0\tCCC\tIII\tNH:i:1\n"
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tGGG\tIII\tNH:i:1\n"
        "reverse-unmapped\t20\t*\t0\t0\t*\t*\t0\t0\tTTT\tIII\tNH:i:1\n"
        "multi\t0\tchr1\t3\t255\t3M\t*\t0\t0\tAAA\tIII\tNH:i:2\n"
    )

    result = subprocess.run(
        ["awk", program, str(alignments)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert [line.split("\t", 1)[0] for line in result.stdout.splitlines()] == [
        "forward",
        "reverse",
    ]

    no_unique = tmp_path / "no-unique.sam"
    no_unique.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tGGG\tIII\tNH:i:1\n"
        "multi\t0\tchr1\t3\t255\t3M\t*\t0\t0\tAAA\tIII\tNH:i:2\n"
    )
    result = subprocess.run(
        ["awk", program, str(no_unique)], capture_output=True, text=True
    )
    assert result.returncode == 1
    assert result.stdout == ""


def test_segemehl_mapping_does_not_pass_a_bare_extension_score_option():
    mapping = rule_body(RULES / "mapping.smk", "map")

    # Segemehl's -e/--extensionscore option requires an integer. The historical
    # bare `-e -d ...` consumed the database flag as that value instead of
    # selecting the documented default extension score.
    assert "segemehl.x -e " not in mapping
    assert (
        "segemehl.x -d {input.genome:q} -i {input.genomeSegemehlIndex:q} "
        "-q {input.fastq:q}"
    ) in mapping


def test_reparation_bam_index_is_adjacent_and_alphanumeric_conditions_render():
    rules = RULES / "reparation.smk"
    adapter = rule_body(rules, "prepareReparationAnnotation")
    rule = rule_body(rules, "reparation")

    assert 'adapter=str(SCRIPTS / "prepare_reparation_annotation.py")' in adapter
    assert 'adapter_deps=[str(SCRIPTS / "gff_utils.py")]' in adapter
    assert (
        "python3 {input.adapter:q} -a {input.annotation:q} -o {output:q}"
        in adapter
    )
    assert "gtf=rules.prepareReparationAnnotation.output" in rule
    assert "gtf=rules.checkAnnotation.output" not in rule
    assert 'bam="maplink/RIBO-{condition}-{replicate}.bam"' in rule
    assert 'bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai"' in rule
    assert "{condition," not in rule
    assert "--output-dir {params.prefix:q}" in rule
    assert "--bai {input.bamindex:q}" in rule


def test_shell_directives_quote_paths_and_track_project_scripts_as_inputs():
    """Keep unsafe interpolation and untracked project code out of shell blocks."""

    unquoted = re.compile(
        r"\{(?:input|output|params|log|wildcards|rules\.)[^{}:]*"
        r"(?:\[[^{}]+\])?\}"
    )
    failures = []
    for path in sorted(RULES.glob("*.smk")):
        for number, block in enumerate(shell_blocks(path.read_text()), start=1):
            matches = [
                match
                for match in unquoted.findall(block)
                # Integer artifact sizes are fixed source constants rather than
                # path or user-controlled interpolation.
                if match != "{params.size}"
            ]
            for untracked_code in ("{SCRIPTS", "{params.script"):
                if untracked_code in block:
                    matches.append(untracked_code)
            if matches:
                failures.append(f"{path.name} shell #{number}: {matches}")

    assert failures == []


def snakemake_command():
    executable = shutil.which("snakemake")
    if executable:
        return [executable]

    conda = shutil.which("conda")
    if conda:
        command = [conda, "run", "-n", "snakemake", "snakemake"]
        probe = subprocess.run(
            [*command, "--version"], capture_output=True, text=True
        )
        if probe.returncode == 0:
            return command
    pytest.skip("Snakemake is unavailable for the shell-rendering dry run")


def test_supported_dag_renders_spaced_paths_and_list_arguments_safely(
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    command = snakemake_command()
    inputs = tmp_path / "external inputs with spaces"
    reads = inputs / "reads with spaces"
    reads.mkdir(parents=True)

    spaced_samples = samples.copy()
    for index, row in spaced_samples.iterrows():
        for column in ("fastqFile", "fastqFile2"):
            source = row[column]
            if not isinstance(source, str) or not source:
                continue
            destination = reads / Path(source).name
            shutil.copyfile(source, destination)
            spaced_samples.at[index, column] = str(destination)

    sample_path = inputs / "sample sheet with spaces.tsv"
    spaced_samples.fillna("").to_csv(sample_path, sep="\t", index=False)
    genome_path = inputs / "reference genome with spaces.fa"
    annotation_path = inputs / "reference annotation with spaces.gff"
    shutil.copyfile(genome_file, genome_path)
    shutil.copyfile(annotation_file, annotation_path)

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome_path),
            "annotation": str(annotation_path),
            "samples": str(sample_path),
        }
    )
    config["workflowSettings"]["stages"] = ["differential_expression"]
    config["predictionSettings"]["deepribo"] = "off"
    config["differentialExpressionSettings"].update(
        {"contrasts": ["B-A"], "features": ["CDS", "custom feature"]}
    )

    config_path = inputs / "workflow config with spaces.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    workdir = tmp_path / "workflow result with spaces"
    workdir.mkdir()
    result = subprocess.run(
        [
            *command,
            "--dry-run",
            "--printshellcmds",
            "--cores",
            "1",
            "--snakefile",
            str(REPO / "workflow" / "Snakefile"),
            "--directory",
            str(workdir),
            "--configfile",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "XDG_CACHE_HOME": str(tmp_path / ".cache")},
    )
    rendered = result.stdout + result.stderr

    assert result.returncode == 0, rendered
    assert shlex.quote(str(spaced_samples.iloc[0]["fastqFile"])) in rendered
    assert "--use_features CDS 'custom feature'" in rendered
    assert rendered.count("touch contrasts/B-A") == 1

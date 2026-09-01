"""Focused regression tests for FastQC rule command construction."""

from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
QC_RULES = (REPO / "workflow" / "rules" / "qc.smk").read_text()


def rule_body(name):
    """Return one rule without relying on Snakemake being installed."""
    start = QC_RULES.index(f"rule {name}:")
    next_rule = QC_RULES.find("\nrule ", start + 1)
    return QC_RULES[start:] if next_rule == -1 else QC_RULES[start:next_rule]


def test_paired_trimmed_fastqc_processes_each_mate_once():
    rule = rule_body("fastqctrimmed_paired")
    commands = [
        line.strip()
        for line in rule.splitlines()
        if line.strip().startswith("fastqc ")
    ]

    assert commands == [
        "fastqc -o {params.outdir:q} -t {threads} {input.reads1:q}",
        "fastqc -o {params.outdir:q} -t {threads} {input.reads2:q}",
    ]
    assert "{input}" not in rule


def test_paired_trimmed_fastqc_publishes_named_outputs_and_logs():
    rule = rule_body("fastqctrimmed_paired")

    assert 'outdir="qc/2trimmed"' in rule
    assert '"logs/{method}-{condition}-{replicate}-trimmed-fastqc.log"' in rule
    assert "exec > {log:q} 2>&1" in rule
    for generated, output in (
        ("params.html1", "output.html1"),
        ("params.zip1", "output.zip1"),
        ("params.html2", "output.html2"),
        ("params.zip2", "output.zip2"),
    ):
        assert f"mv {{{generated}:q}} {{{output}:q}}" in rule

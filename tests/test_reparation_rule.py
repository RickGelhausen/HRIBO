"""Regression tests for the REPARATION execution and artifact contract."""

import importlib.util
import os
import signal
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parent.parent
RULE = (REPO / "workflow" / "rules" / "reparation.smk").read_text()
RUNNER = REPO / "workflow" / "scripts" / "run_reparation.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("run_reparation_under_test", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def reparation_run(tmp_path):
    """Create spaced input/output paths and a controllable fake engine."""
    tool_dir = tmp_path / "fake tools"
    tool_dir.mkdir()
    engine = tool_dir / "reparation fake.py"
    engine.write_text(
        r'''#!/usr/bin/env python3
import base64
import os
import re
import shutil
import sys
import time
from pathlib import Path


def option(name):
    return Path(sys.argv[sys.argv.index(name) + 1])


mode = os.environ.get("FAKE_REPARATION_MODE", "valid_png")
genome = option("-g")
gtf = option("-gtf")
database = option("-db")
bam = option("-bam")
work_dir = option("-wdir")

for path in (genome, gtf, database, bam, work_dir):
    assert re.fullmatch(r"[A-Za-z0-9_./-]+", str(path)), path
assert genome.name == "genome.fa" and genome.is_symlink()
assert gtf.name == "annotation.gtf" and gtf.is_symlink()
assert bam.name == "reads.bam" and bam.is_symlink()
assert Path(str(bam) + ".bai").is_symlink()
assert database.name == "protein_db.fasta" and database.is_symlink()

Path(str(database) + ".pin").write_text("isolated BLAST sidecar\n")
if mode == "pause_in_staging":
    (work_dir.parent / "large-staging-artifact.bin").write_bytes(b"x" * 1024 * 1024)
    notification = Path(os.environ["FAKE_REPARATION_NOTIFICATION"])
    with notification.open("w") as handle:
        handle.write(str(os.getpid()))
        handle.flush()
        os.fsync(handle.fileno())
    while True:
        time.sleep(1)
if mode == "fail":
    print("simulated REPARATION failure", file=sys.stderr)
    sys.exit(23)

if work_dir.exists():
    shutil.rmtree(str(work_dir))
work_dir.mkdir(parents=True)

header = (
    "ORF_locus\tstrand\tlength\tstart_codon\tribo_count\tribo_rpkm\t"
    "ribo_coverage\tSD_score\tSD_pos\tprob\tORF_type\tReference\t"
    "Distance_from_aTIS\n"
)
contig, start, stop, strand = "chr1", 10, 39, "+"
if mode == "unknown_contig":
    contig = "missing"
elif mode == "plus_out_of_bounds":
    start, stop = 70, 99
elif mode == "plus_exact_boundary":
    start, stop = 68, 97
elif mode == "minus_out_of_bounds":
    start, stop, strand = 3, 32, "-"
elif mode == "minus_exact_boundary":
    start, stop, strand = 4, 33, "-"
elif mode == "colon_contig":
    contig = "chr:alpha"
locus = f"{contig}:{start}-{stop}"
row = (
    f"{locus}\t{strand}\t30\tATG\t12\t3.5\t0.8\t-2\t4\t0.95\t"
    "Annotated\tgeneA\tNA\n"
)
if mode == "bad_table":
    row = "chr1:10-39\t+\t30\tATG\t12\t3.5\t1.8\t-2\t4\t0.95\tAnnotated\tgeneA\tNA\n"
zero_predictions = mode == "zero"
table_rows = "" if zero_predictions else row
if mode == "duplicate_table":
    table_rows += row
(work_dir / "Predicted_ORFs.txt").write_text(header + table_rows)

bed_header = 'track type=bed name="Predicted_ORFs" description="" visibility=1\n'
if strand == "+":
    bed_start, bed_stop, thick_start, thick_stop = start - 1, stop + 3, start - 1, stop
else:
    bed_start, bed_stop, thick_start, thick_stop = start - 4, stop, start - 1, stop
bed_row = (
    f"{contig}\t{bed_start}\t{bed_stop}\t{locus}\t1\t{strand}\t"
    f"{thick_start}\t{thick_stop}\t7,7,255\n"
)
if mode == "bed_count":
    bed_row = ""
elif mode == "bed_coordinates":
    bed_row = "chr1\t8\t42\tchr1:10-39\t1\t+\t9\t39\t7,7,255\n"
(work_dir / "Predicted_ORFs.bed").write_text(
    bed_header + ("" if zero_predictions else bed_row)
)

fasta = f">generic|{locus}|start codon:ATG strand:{strand} length:30\nMABCDEFGHI\n"
if zero_predictions or mode == "fasta_count":
    fasta = ""
elif mode == "fasta_metadata":
    fasta = ">generic|chr1:10-39|start codon:GTG strand:+ length:30\nMABCDEFGHI\n"
elif mode == "fasta_sequence_length":
    fasta = ">generic|chr1:10-39|start codon:ATG strand:+ length:30\nMABCDEFGH\n"
(work_dir / "Predicted_ORFs.fasta").write_text(fasta)

def minimal_pdf():
    objects = (
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 1 1] /Contents 4 0 R >>",
        b"<< /Length 0 >>\nstream\n\nendstream",
    )
    content = b"%PDF-1.4\n"
    offsets = [0]
    for number, body in enumerate(objects, 1):
        offsets.append(len(content))
        content += str(number).encode() + b" 0 obj\n" + body + b"\nendobj\n"
    xref = len(content)
    content += b"xref\n0 5\n0000000000 65535 f \n"
    for offset in offsets[1:]:
        content += ("%010d 00000 n \n" % offset).encode()
    content += (
        b"trailer\n<< /Size 5 /Root 1 0 R >>\nstartxref\n"
        + str(xref).encode()
        + b"\n%%EOF\n"
    )
    return content


pdf = minimal_pdf()
for name in (
    "metagene_profile.pdf",
    "PR_and_ROC_curve.pdf",
    "variable_importance.pdf",
    "S_Curve.pdf",
):
    (work_dir / name).write_bytes(pdf)
if mode == "bad_pdf":
    (work_dir / "S_Curve.pdf").write_text("not a PDF\n")
if mode == "missing":
    (work_dir / "variable_importance.pdf").unlink()
if mode == "result_symlink":
    (work_dir / "unexpected-link").symlink_to(database)
if mode == "reserved_receipt":
    (work_dir / ".complete").write_text("engine-controlled receipt\n")

offsets = "## generated by psite\n\nlength\tp_offset\n28\t13\ndefault\t13\n"
if mode == "bad_offsets":
    offsets = "length\tp_offset\n"
elif mode == "offset_at_length":
    offsets = "length\tp_offset\n28\t28\n"
(work_dir / "p_site_offsets.txt").write_text(offsets)

png = work_dir / "p_site_offset.png"
if mode == "empty_png":
    png.touch()
elif mode == "bad_png":
    png.write_text("not a PNG\n")
elif mode != "no_png":
    png.write_bytes(
        base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
        )
    )
'''
    )
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR)

    input_dir = tmp_path / "input files"
    input_dir.mkdir()
    paths = {}
    for key, name in (
        ("genome", "reference genome.fa"),
        ("gtf", "reference annotation.gtf"),
        ("database", "protein database.fasta"),
        ("bam", "mapped reads.bam"),
        ("bai", "mapped reads.bam.bai"),
    ):
        paths[key] = input_dir / name
        paths[key].write_text("fixture\n")
    paths["genome"].write_text(
        ">chr1 reference description\n"
        + "A" * 100
        + "\n>chr:alpha colon-bearing identifier\n"
        + "C" * 100
        + "\n"
    )

    workflow_root = tmp_path / "workflow root"
    output_dir = workflow_root / "reparation" / "published result"
    output_dir.parent.mkdir(parents=True)
    command = [
        sys.executable,
        str(RUNNER),
        "--engine",
        str(engine),
        "--genome",
        str(paths["genome"]),
        "--gtf",
        str(paths["gtf"]),
        "--database",
        str(paths["database"]),
        "--bam",
        str(paths["bam"]),
        "--bai",
        str(paths["bai"]),
        "--output-dir",
        str(output_dir),
        "--threads",
        "3",
    ]
    return command, paths, output_dir


def execute(fixture, mode):
    command, _, _ = fixture
    output_dir = fixture[2]
    environment = {**os.environ, "FAKE_REPARATION_MODE": mode}
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=environment,
        cwd=output_dir.parent.parent,
    )


def transaction_artifact(output_dir, suffix):
    return output_dir.parent / f".{output_dir.name}_{suffix}"


def result_snapshot(output_dir):
    return {
        path.relative_to(output_dir): path.read_bytes()
        for path in output_dir.rglob("*")
        if path.is_file() and path.name != ".complete"
    }


def production_rule_run(reparation_run, snakemake_command):
    """Stage the real Reparation rule around the fixture's fake engine."""
    command, paths, fixture_output = reparation_run
    workflow_root = fixture_output.parent.parent
    output_dir = workflow_root / "reparation/A-1"

    fake_bin = workflow_root / "production-fake-bin"
    fake_bin.mkdir()
    (fake_bin / "reparation.pl").symlink_to(command[command.index("--engine") + 1])

    staged_inputs = {
        workflow_root / "genomes/genome.fa": paths["genome"],
        workflow_root / "annotation/annotation.gff": paths["gtf"],
        workflow_root / "reparation/annotation.gtf": paths["gtf"],
        workflow_root / "uniprotDB/uniprot_sprot.fasta": paths["database"],
        workflow_root / "maplink/RIBO-A-1.bam": paths["bam"],
        workflow_root / "maplink/RIBO-A-1.bam.bai": paths["bai"],
    }
    for destination, source in staged_inputs.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source.read_bytes())

    snakefile = workflow_root / "ProductionReparationSnakefile"
    snakefile.write_text(
        f'''from pathlib import Path
import pandas as pd

SCRIPTS = Path({str(REPO / "workflow/scripts")!r})
samples = pd.DataFrame(
    [{{"method": "RIBO", "condition": "A", "replicate": "1"}}]
)


rule retrieveGenome:
    output:
        "genomes/genome.fa"


rule checkAnnotation:
    output:
        "annotation/annotation.gff"


include: {str(REPO / "workflow/rules/reparation.smk")!r}
'''
    )
    invocation = [
        *snakemake_command,
        (output_dir / ".complete").relative_to(workflow_root).as_posix(),
        "--snakefile",
        str(snakefile),
        "--cores",
        "1",
        "--resources",
        "reparation_instances=1",
        "--allowed-rules",
        "reparation",
        "--printshellcmds",
        "--show-failed-logs",
    ]
    environment = {
        **os.environ,
        "PATH": os.pathsep.join([str(fake_bin), os.environ.get("PATH", "")]),
        "XDG_CACHE_HOME": str(workflow_root / ".cache"),
    }
    return invocation, environment, workflow_root, output_dir


def test_reparation_rule_delegates_to_the_checked_runner():
    assert 'runner=workflow.source_path("../scripts/run_reparation.py")' in RULE
    assert "python3 {input.runner:q}" in RULE
    assert 'receipt=ensure(' in RULE
    assert '"reparation/{condition}-{replicate}/.complete", non_empty=True' in RULE
    assert "receipt=rules.reparation.output.receipt" in RULE
    assert 'os.path.dirname(input.receipt), "Predicted_ORFs.txt"' in RULE
    assert "-i {params.orfs:q}" in RULE
    assert "p_site_offset.png" not in RULE
    for argument in (
        "--genome {input.genome:q}",
        "--gtf {input.gtf:q}",
        "--database {input.db:q}",
        "--bam {input.bam:q}",
        "--bai {input.bamindex:q}",
        "--output-dir {params.prefix:q}",
        "--threads {threads:q}",
    ):
        assert argument in RULE
    assert "|| true" not in RULE
    assert "touch {output" not in RULE
    assert "update(" not in RULE
    assert "directory(" not in RULE


def test_reparation_stages_safe_aliases_and_replaces_stale_results(reparation_run):
    _, paths, output_dir = reparation_run
    output_dir.mkdir()
    stale = output_dir / "stale artifact.txt"
    stale.write_text("old result\n")

    result = execute(reparation_run, "valid_png")

    assert result.returncode == 0, result.stderr
    assert not stale.exists()
    assert (output_dir / "Predicted_ORFs.txt").is_file()
    assert (output_dir / "p_site_offset.png").stat().st_size > 0
    assert (output_dir / ".complete").read_bytes() == (
        b"HRIBO REPARATION publication v1\n"
    )
    assert not Path(str(paths["database"]) + ".pin").exists()
    for suffix in (
        "candidate",
        "backup",
        "transaction",
        "transaction.tmp",
        "staging",
        "staging.tmp",
    ):
        assert not transaction_artifact(output_dir, suffix).exists()


def test_reparation_reclaims_only_its_own_stage_after_process_group_sigkill(
    reparation_run, tmp_path
):
    command, _, output_dir = reparation_run
    runner = load_runner()
    first_temporary_base = tmp_path / "safe_tmp_one"
    second_temporary_base = tmp_path / "safe_tmp_two"
    first_temporary_base.mkdir(mode=0o700)
    second_temporary_base.mkdir(mode=0o700)
    assert runner.SAFE_PATH.fullmatch(str(first_temporary_base.resolve()))
    assert runner.SAFE_PATH.fullmatch(str(second_temporary_base.resolve()))

    notification = tmp_path / "engine-in-staging"
    environment = {
        **os.environ,
        "TMPDIR": str(first_temporary_base),
        "FAKE_REPARATION_MODE": "pause_in_staging",
        "FAKE_REPARATION_NOTIFICATION": str(notification),
    }
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=output_dir.parent.parent,
        env=environment,
        start_new_session=True,
    )
    try:
        deadline = time.monotonic() + 30
        while not notification.exists() and process.poll() is None:
            if time.monotonic() >= deadline:
                break
            time.sleep(0.01)
        assert notification.exists(), "fake engine never populated its staging tree"
        stage_root = runner._staging_path(
            output_dir, str(first_temporary_base.resolve())
        )
        assert (stage_root / "large-staging-artifact.bin").stat().st_size == 1024 * 1024
        assert transaction_artifact(output_dir, "staging").is_file()

        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate(timeout=30)
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=10)

    assert process.returncode != 0, stdout + stderr
    assert stage_root.is_dir()

    # A sibling with another output's digest represents a concurrent stage in
    # the same private namespace and must never be swept by this recovery.
    other_output = output_dir.parent / "another sample"
    other_stage = runner._staging_path(
        other_output, str(first_temporary_base.resolve())
    )
    other_stage.mkdir(mode=0o700)
    (other_stage / runner.STAGING_MARKER).write_bytes(
        runner._staging_marker_content(other_output)
    )
    other_sentinel = other_stage / "must-remain.bin"
    other_sentinel.write_bytes(b"unrelated staging data\n")

    restarted = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
        env={
            **os.environ,
            "TMPDIR": str(second_temporary_base),
            "FAKE_REPARATION_MODE": "valid_png",
        },
    )

    assert restarted.returncode == 0, restarted.stderr
    assert not stage_root.exists()
    assert other_sentinel.read_bytes() == b"unrelated staging data\n"
    assert not transaction_artifact(output_dir, "staging").exists()
    assert not transaction_artifact(output_dir, "staging.tmp").exists()
    assert (output_dir / ".complete").is_file()


def test_reparation_refuses_a_symlinked_stale_staging_tree(
    reparation_run, tmp_path
):
    command, _, output_dir = reparation_run
    runner = load_runner()
    temporary_base = tmp_path / "safe_tmp"
    temporary_base.mkdir(mode=0o700)
    namespace = runner._staging_namespace(str(temporary_base.resolve()))
    namespace.mkdir(mode=0o700)
    stage_root = runner._staging_path(output_dir, str(temporary_base.resolve()))
    unrelated = tmp_path / "unrelated-staging-target"
    unrelated.mkdir()
    sentinel = unrelated / "must remain.txt"
    sentinel.write_text("untouched\n")
    stage_root.symlink_to(unrelated, target_is_directory=True)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
        env={
            **os.environ,
            "TMPDIR": str(temporary_base),
            "FAKE_REPARATION_MODE": "valid_png",
        },
    )

    assert result.returncode == 1
    assert "staging directory is not a real directory" in result.stderr
    assert stage_root.is_symlink()
    assert sentinel.read_text() == "untouched\n"
    assert not output_dir.exists()


def test_reparation_refuses_a_permissive_staging_locator(reparation_run, tmp_path):
    command, _, output_dir = reparation_run
    runner = load_runner()
    temporary_base = tmp_path / "safe_tmp"
    temporary_base.mkdir(mode=0o700)
    namespace = runner._staging_namespace(str(temporary_base.resolve()))
    namespace.mkdir(mode=0o700)
    stage_root = runner._staging_path(output_dir, str(temporary_base.resolve()))
    stage_root.mkdir(mode=0o700)
    marker = stage_root / runner.STAGING_MARKER
    marker.write_bytes(runner._staging_marker_content(output_dir))
    marker.chmod(0o600)
    sentinel = stage_root / "must-remain.bin"
    sentinel.write_bytes(b"not authorized by a private locator\n")
    locator = transaction_artifact(output_dir, "staging")
    locator.write_bytes(runner._staging_locator_content(output_dir, stage_root))
    locator.chmod(0o644)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
        env={
            **os.environ,
            "TMPDIR": str(temporary_base),
            "FAKE_REPARATION_MODE": "valid_png",
        },
    )

    assert result.returncode == 1
    assert "staging locator must be owned" in result.stderr
    assert "mode 0600" in result.stderr
    assert sentinel.read_bytes() == b"not authorized by a private locator\n"
    assert locator.is_file()
    assert not output_dir.exists()


def test_reparation_refuses_a_locator_below_an_unsafe_temporary_base(
    reparation_run, tmp_path
):
    command, _, output_dir = reparation_run
    runner = load_runner()
    selected_base = tmp_path / "selected_safe_tmp"
    selected_base.mkdir(mode=0o700)
    unsafe_base = tmp_path / "unsafe_tmp"
    unsafe_base.mkdir(mode=0o700)
    unsafe_base.chmod(0o777)
    namespace = runner._staging_namespace(str(unsafe_base.resolve()))
    namespace.mkdir(mode=0o700)
    stage_root = runner._staging_path(output_dir, str(unsafe_base.resolve()))
    stage_root.mkdir(mode=0o700)
    marker = stage_root / runner.STAGING_MARKER
    marker.write_bytes(runner._staging_marker_content(output_dir))
    marker.chmod(0o600)
    sentinel = stage_root / "must-remain.bin"
    sentinel.write_bytes(b"unsafe-parent staging data\n")
    locator = transaction_artifact(output_dir, "staging")
    locator.write_bytes(runner._staging_locator_content(output_dir, stage_root))
    locator.chmod(0o600)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
        env={
            **os.environ,
            "TMPDIR": str(selected_base),
            "FAKE_REPARATION_MODE": "valid_png",
        },
    )

    assert result.returncode == 1
    assert "temporary base is writable and not sticky" in result.stderr
    assert sentinel.read_bytes() == b"unsafe-parent staging data\n"
    assert locator.is_file()
    assert not output_dir.exists()


def test_reparation_propagates_exit_and_preserves_stale_results(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    stale = output_dir / "stale artifact.txt"
    stale.write_text("old result\n")
    receipt = output_dir / ".complete"
    receipt.write_text("old completion receipt\n")
    temporary_receipt = output_dir / ".complete.tmp"
    temporary_receipt.write_text("interrupted receipt publication\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert "simulated REPARATION failure" in result.stderr
    assert stale.read_text() == "old result\n"
    assert not receipt.exists()
    assert not temporary_receipt.exists()
    assert sorted(path.name for path in output_dir.iterdir()) == [stale.name]


def test_reparation_recovers_an_original_moved_before_shutdown(reparation_run):
    _, _, output_dir = reparation_run
    backup = transaction_artifact(output_dir, "backup")
    backup.mkdir()
    (backup / "old.txt").write_text("old result\n")
    candidate = transaction_artifact(output_dir, "candidate")
    candidate.mkdir()
    (candidate / "partial.txt").write_text("uncommitted result\n")
    transaction_artifact(output_dir, "transaction").write_text("present\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert (output_dir / "old.txt").read_text() == "old result\n"
    assert not backup.exists()
    assert not candidate.exists()
    assert not transaction_artifact(output_dir, "transaction").exists()


def test_reparation_rolls_back_a_new_tree_installed_before_commit(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    (output_dir / "new.txt").write_text("uncommitted result\n")
    backup = transaction_artifact(output_dir, "backup")
    backup.mkdir()
    (backup / "old.txt").write_text("old result\n")
    transaction_artifact(output_dir, "transaction").write_text("present\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert (output_dir / "old.txt").read_text() == "old result\n"
    assert not (output_dir / "new.txt").exists()
    assert not backup.exists()
    assert not transaction_artifact(output_dir, "transaction").exists()


def test_reparation_recovers_a_shutdown_before_the_original_move(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    (output_dir / "old.txt").write_text("old result\n")
    candidate = transaction_artifact(output_dir, "candidate")
    candidate.mkdir()
    (candidate / "complete.txt").write_text("uncommitted result\n")
    transaction_artifact(output_dir, "transaction").write_text("present\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert (output_dir / "old.txt").read_text() == "old result\n"
    assert not candidate.exists()
    assert not transaction_artifact(output_dir, "transaction").exists()


@pytest.mark.parametrize("interruption", ["before_install", "after_install"])
def test_reparation_rolls_back_an_interrupted_first_publication(
    reparation_run, interruption
):
    _, _, output_dir = reparation_run
    candidate = transaction_artifact(output_dir, "candidate")
    if interruption == "before_install":
        candidate.mkdir()
        (candidate / "partial.txt").write_text("uncommitted result\n")
    else:
        output_dir.mkdir()
        (output_dir / "partial.txt").write_text("uncommitted result\n")
    transaction_artifact(output_dir, "transaction").write_text("absent\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert not output_dir.exists()
    assert not candidate.exists()
    assert not transaction_artifact(output_dir, "transaction").exists()


def test_reparation_treats_post_commit_cleanup_failure_as_success(
    tmp_path, monkeypatch, capsys
):
    runner = load_runner()
    staged = tmp_path / "staged"
    staged.mkdir()
    (staged / "new.txt").write_text("new result\n")
    output_dir = tmp_path / "reparation" / "sample"
    output_dir.mkdir(parents=True)
    (output_dir / "old.txt").write_text("old result\n")
    monkeypatch.setattr(runner, "_validate_result", lambda *_: None)
    monkeypatch.setattr(runner, "_sync_tree", lambda *_: None)
    real_sync = runner._sync_directory
    calls = 0

    def fail_after_commit(path):
        nonlocal calls
        calls += 1
        # Candidate prep, marker creation, old move, and new move have already
        # been synchronized. The fifth call follows marker removal.
        if calls == 5:
            raise OSError("simulated post-commit fsync failure")
        real_sync(path)

    monkeypatch.setattr(runner, "_sync_directory", fail_after_commit)

    runner._publish(staged, output_dir, {})

    assert (output_dir / "new.txt").read_text() == "new result\n"
    backup = transaction_artifact(output_dir, "backup")
    assert (backup / "old.txt").read_text() == "old result\n"
    assert not transaction_artifact(output_dir, "transaction").exists()
    assert "output was committed" in capsys.readouterr().err

    monkeypatch.setattr(runner, "_sync_directory", real_sync)
    runner._recover_publication(output_dir)
    assert (output_dir / "new.txt").read_text() == "new result\n"
    assert not backup.exists()


def test_reparation_finishes_cleanup_after_a_committed_publication(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    (output_dir / "new.txt").write_text("new result\n")
    backup = transaction_artifact(output_dir, "backup")
    backup.mkdir()
    (backup / "old.txt").write_text("old result\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert (output_dir / "new.txt").read_text() == "new result\n"
    assert not backup.exists()


def test_reparation_recovers_a_single_legacy_random_backup(reparation_run):
    _, _, output_dir = reparation_run
    backup = output_dir.parent / f".{output_dir.name}_backup_deadbeef"
    backup.mkdir()
    (backup / "old.txt").write_text("old result\n")
    publish_root = output_dir.parent / f".{output_dir.name}_publish_deadbeef"
    (publish_root / "result").mkdir(parents=True)

    result = execute(reparation_run, "fail")

    assert result.returncode == 23
    assert (output_dir / "old.txt").read_text() == "old result\n"
    assert not backup.exists()
    assert not publish_root.exists()


def test_reparation_refuses_ambiguous_legacy_backups(reparation_run):
    _, _, output_dir = reparation_run
    backups = [
        output_dir.parent / f".{output_dir.name}_backup_{suffix}"
        for suffix in ("first", "second")
    ]
    for backup in backups:
        backup.mkdir()
        (backup / "old.txt").write_text(f"{backup.name}\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 1
    assert "multiple legacy publication backups" in result.stderr
    assert not output_dir.exists()
    assert all(backup.exists() for backup in backups)


def test_reparation_refuses_a_symlinked_transaction_backup(reparation_run):
    _, _, output_dir = reparation_run
    target = output_dir.parent / "unrelated backup target"
    target.mkdir()
    sentinel = target / "keep.txt"
    sentinel.write_text("unchanged\n")
    transaction_artifact(output_dir, "backup").symlink_to(
        target, target_is_directory=True
    )
    transaction_artifact(output_dir, "transaction").write_text("present\n")

    result = execute(reparation_run, "fail")

    assert result.returncode == 1
    assert "publication backup is not a real directory" in result.stderr
    assert sentinel.read_text() == "unchanged\n"
    assert transaction_artifact(output_dir, "backup").is_symlink()


def test_reparation_refuses_a_symlinked_completion_receipt(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    stale = output_dir / "old.txt"
    stale.write_text("old result\n")
    unrelated = output_dir.parent / "unrelated receipt target"
    unrelated.write_text("must remain unchanged\n")
    receipt = output_dir / ".complete"
    receipt.symlink_to(unrelated)

    result = execute(reparation_run, "fail")

    assert result.returncode == 1
    assert "completion receipt is not a regular file" in result.stderr
    assert receipt.is_symlink()
    assert unrelated.read_text() == "must remain unchanged\n"
    assert stale.read_text() == "old result\n"


def test_reparation_rejects_an_output_directory_symlink(reparation_run):
    command, _, output_dir = reparation_run
    target = output_dir.parent / "symlink target"
    target.mkdir()
    sentinel = target / "must remain.txt"
    sentinel.write_text("unchanged\n")
    output_dir.symlink_to(target, target_is_directory=True)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
    )

    assert result.returncode != 0
    assert "not a file or symlink" in result.stderr
    assert output_dir.is_symlink()
    assert sentinel.read_text() == "unchanged\n"


def test_reparation_rejects_a_non_directory_output_target(reparation_run):
    _, _, output_dir = reparation_run
    output_dir.write_text("must remain a file\n")

    result = execute(reparation_run, "valid_png")

    assert result.returncode != 0
    assert "not a file or symlink" in result.stderr
    assert output_dir.read_text() == "must remain a file\n"


def test_reparation_refuses_the_filesystem_root(reparation_run):
    command, _, output_dir = reparation_run
    command[command.index("--output-dir") + 1] = "/"

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
    )

    assert result.returncode != 0
    assert "refusing to use the filesystem root" in result.stderr


def test_reparation_refuses_the_workflow_root(reparation_run):
    command, _, output_dir = reparation_run
    workflow_root = output_dir.parent.parent
    command[command.index("--output-dir") + 1] = "."

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=workflow_root,
    )

    assert result.returncode != 0
    assert "must be a direct child" in result.stderr


def test_reparation_refuses_a_symlinked_output_parent(reparation_run):
    command, _, output_dir = reparation_run
    workflow_root = output_dir.parent.parent
    output_dir.parent.rmdir()
    target = workflow_root / "real result parent"
    target.mkdir()
    output_dir.parent.symlink_to(target, target_is_directory=True)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=workflow_root,
    )

    assert result.returncode != 0
    assert "symlinked reparation directory" in result.stderr
    assert output_dir.parent.is_symlink()


def test_reparation_refuses_an_output_containing_an_input(reparation_run):
    command, _, output_dir = reparation_run
    output_dir.mkdir()
    nested_genome = output_dir / "genome.fa"
    nested_genome.write_text("must remain\n")
    command[command.index("--genome") + 1] = str(nested_genome)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        cwd=output_dir.parent.parent,
    )

    assert result.returncode != 0
    assert "contains the genome input" in result.stderr
    assert nested_genome.read_text() == "must remain\n"


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("missing", "required PDF"),
        ("bad_table", "ribo_coverage above 1"),
        ("bed_count", "has 0 records; expected 1"),
        ("bed_coordinates", "coordinates inconsistent with the ORF table"),
        ("duplicate_table", "duplicate ORF"),
        ("fasta_count", "empty protein FASTA"),
        ("fasta_metadata", "codon or length inconsistent with the ORF table"),
        ("fasta_sequence_length", "codon or length inconsistent with the ORF table"),
        ("bad_pdf", "S_Curve.pdf is not a valid PDF"),
        ("bad_offsets", "has no offset data rows"),
        ("offset_at_length", "offset outside the read length"),
        ("bad_png", "is not a valid PNG"),
        ("result_symlink", "symlink or special file"),
        ("reserved_receipt", "reserved publication file"),
    ],
)
def test_reparation_rejects_missing_or_malformed_artifacts(
    reparation_run, mode, message
):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    stale = output_dir / "stale artifact.txt"
    stale.write_text("old result\n")

    result = execute(reparation_run, mode)

    assert result.returncode != 0
    assert message in result.stderr
    assert stale.read_text() == "old result\n"
    assert sorted(path.name for path in output_dir.iterdir()) == [stale.name]


def test_reparation_accepts_a_valid_zero_prediction_result(reparation_run):
    result = execute(reparation_run, "zero")
    _, _, output_dir = reparation_run

    assert result.returncode == 0, result.stderr
    assert (output_dir / "Predicted_ORFs.txt").read_text().count("\n") == 1
    assert (output_dir / "Predicted_ORFs.bed").read_text().count("\n") == 1
    assert (output_dir / "Predicted_ORFs.fasta").stat().st_size == 0


@pytest.mark.parametrize(
    "mode", ["plus_exact_boundary", "minus_exact_boundary", "colon_contig"]
)
def test_reparation_accepts_stop_inclusive_reference_boundaries(
    reparation_run, mode
):
    result = execute(reparation_run, mode)
    _, _, output_dir = reparation_run

    assert result.returncode == 0, result.stderr
    assert (output_dir / "Predicted_ORFs.txt").is_file()


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("unknown_contig", "unknown genome contig 'missing'"),
        ("plus_out_of_bounds", "stop-inclusive interval 70-102 outside"),
        ("minus_out_of_bounds", "stop-inclusive interval 0-32 outside"),
    ],
)
def test_reparation_rejects_predictions_outside_the_reference(
    reparation_run, mode, message
):
    _, _, output_dir = reparation_run
    output_dir.mkdir()
    sentinel = output_dir / "old.txt"
    sentinel.write_text("old result\n")

    result = execute(reparation_run, mode)

    assert result.returncode == 1
    assert message in result.stderr
    assert sentinel.read_text() == "old result\n"


def test_reparation_rejects_duplicate_reference_identifiers(reparation_run):
    _, paths, output_dir = reparation_run
    paths["genome"].write_text(">chr1 first\nAAAA\n>chr1 second\nCCCC\n")

    result = execute(reparation_run, "valid_png")

    assert result.returncode == 1
    assert "duplicate record 'chr1'" in result.stderr
    assert not output_dir.exists()


def test_reparation_rejects_non_iupac_reference_sequence(reparation_run):
    _, paths, output_dir = reparation_run
    paths["genome"].write_text(">chr1\nAAAA10 bases\n")

    result = execute(reparation_run, "valid_png")

    assert result.returncode == 1
    assert "non-IUPAC sequence data" in result.stderr
    assert not output_dir.exists()


def test_reparation_removes_the_empty_upstream_png_fallback(reparation_run):
    result = execute(reparation_run, "empty_png")
    _, _, output_dir = reparation_run

    assert result.returncode == 0, result.stderr
    assert "empty p_site_offset.png fallback" in result.stderr
    assert not (output_dir / "p_site_offset.png").exists()


def test_reparation_keeps_a_genuine_optional_png(reparation_run):
    result = execute(reparation_run, "valid_png")
    _, _, output_dir = reparation_run

    assert result.returncode == 0, result.stderr
    assert (output_dir / "p_site_offset.png").read_bytes().startswith(b"\x89PNG")


def test_snakemake_preserves_the_result_tree_after_engine_failure(
    reparation_run, snakemake_command
):
    invocation, run_environment, workflow_root, output_dir = production_rule_run(
        reparation_run, snakemake_command
    )

    successful = subprocess.run(
        invocation,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )
    assert successful.returncode == 0, successful.stderr
    optional_png = output_dir / "p_site_offset.png"
    optional_png.write_bytes(b"old optional image\n")
    old_extra = output_dir / "old-only.txt"
    old_extra.write_text("whole old snapshot\n")
    previous = result_snapshot(output_dir)

    failed = subprocess.run(
        [*invocation, "--forcerun", "reparation"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "fail"},
    )

    assert failed.returncode != 0
    assert "simulated REPARATION failure" in failed.stderr
    assert result_snapshot(output_dir) == previous
    assert not (output_dir / ".complete").exists()


def test_missing_required_sibling_fails_downstream_until_receipt_is_invalidated(
    reparation_run, snakemake_command
):
    invocation, run_environment, workflow_root, output_dir = production_rule_run(
        reparation_run, snakemake_command
    )
    successful = subprocess.run(
        invocation,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )
    assert successful.returncode == 0, successful.stderr

    receipt = output_dir / ".complete"
    orfs = output_dir / "Predicted_ORFs.txt"
    orfs.unlink()
    assert receipt.is_file()

    downstream = list(invocation)
    downstream[len(snakemake_command)] = "reparation/A-1.reparation.gff"
    allowed = downstream.index("--allowed-rules")
    del downstream[allowed : allowed + 2]
    failed = subprocess.run(
        downstream,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )

    assert failed.returncode != 0
    assert "Predicted_ORFs.txt" in failed.stdout + failed.stderr
    assert receipt.is_file()
    assert not (workflow_root / "reparation/A-1.reparation.gff").exists()

    # The receipt is the unit Snakemake tracks. Invalidating it schedules the
    # owner, which atomically replaces the incomplete tree and republishes it.
    receipt.unlink()
    repaired = subprocess.run(
        invocation,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )
    assert repaired.returncode == 0, repaired.stderr
    assert orfs.is_file()
    assert receipt.read_bytes() == b"HRIBO REPARATION publication v1\n"


def test_production_rule_recovers_after_whole_process_group_postcommit_sigkill(
    reparation_run, snakemake_command
):
    invocation, run_environment, workflow_root, output_dir = production_rule_run(
        reparation_run, snakemake_command
    )
    successful = subprocess.run(
        invocation,
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )
    assert successful.returncode == 0, successful.stderr

    orfs = output_dir / "Predicted_ORFs.txt"
    orfs.write_text(orfs.read_text() + "old declared snapshot\n")
    optional_png = output_dir / "p_site_offset.png"
    optional_png.write_bytes(b"old optional image\n")
    old_extra = output_dir / "old-only.txt"
    old_extra.write_text("whole old snapshot\n")
    old_snapshot = result_snapshot(output_dir)

    hook_dir = workflow_root / "commit-pause-hook"
    hook_dir.mkdir()
    notification = workflow_root / "runner-committed.pid"
    (hook_dir / "sitecustomize.py").write_text(
        r'''import os
import sys
import time


if (
    os.environ.get("HRIBO_TEST_COMMIT_NOTIFY")
    and sys.argv
    and sys.argv[0].endswith("run_reparation.py")
):
    real_unlink = os.unlink

    def pause_after_transaction_unlink(path, *args, **kwargs):
        result = real_unlink(path, *args, **kwargs)
        if os.path.basename(os.fspath(path)).endswith("_transaction"):
            with open(os.environ["HRIBO_TEST_COMMIT_NOTIFY"], "w") as handle:
                handle.write(str(os.getpid()))
                handle.flush()
                os.fsync(handle.fileno())
            while True:
                time.sleep(1)
        return result

    os.unlink = pause_after_transaction_unlink
'''
    )
    kill_environment = {
        **run_environment,
        "FAKE_REPARATION_MODE": "valid_png",
        "HRIBO_TEST_COMMIT_NOTIFY": str(notification),
        "PYTHONPATH": os.pathsep.join(
            [str(hook_dir), run_environment.get("PYTHONPATH", "")]
        ).rstrip(os.pathsep),
    }
    process = subprocess.Popen(
        [*invocation, "--forcerun", "reparation"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=workflow_root,
        env=kill_environment,
        start_new_session=True,
    )
    try:
        deadline = time.monotonic() + 60
        while not notification.exists() and process.poll() is None:
            if time.monotonic() >= deadline:
                break
            time.sleep(0.01)
        assert notification.exists(), "runner never reached its publication commit"
        assert int(notification.read_text()) != process.pid
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate(timeout=120)
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=10)

    assert process.returncode != 0, stdout + stderr
    assert not (output_dir / ".complete").exists()
    assert not old_extra.exists()
    assert b"old declared snapshot" not in orfs.read_bytes()
    assert optional_png.read_bytes().startswith(b"\x89PNG")
    committed_snapshot = result_snapshot(output_dir)
    assert committed_snapshot != old_snapshot

    # Unlike update(directory(...)), a killed scheduler leaves no private
    # Snakemake backup that can make jobs.prepare() fail before the runner gets
    # a chance to reconcile its own transaction.
    assert not (workflow_root / ".snakemake/backups").exists()
    backup = transaction_artifact(output_dir, "backup")
    assert backup.is_dir()
    unlocked = subprocess.run(
        [*invocation, "--unlock"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env=run_environment,
    )
    assert unlocked.returncode == 0, unlocked.stderr

    recovered = subprocess.run(
        [*invocation, "--rerun-incomplete", "--forcerun", "reparation"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "fail"},
    )
    assert recovered.returncode != 0
    assert "simulated REPARATION failure" in recovered.stderr
    assert result_snapshot(output_dir) == committed_snapshot
    assert not (output_dir / ".complete").exists()
    assert not backup.exists()

    restarted = subprocess.run(
        [*invocation, "--rerun-incomplete"],
        capture_output=True,
        text=True,
        cwd=workflow_root,
        env={**run_environment, "FAKE_REPARATION_MODE": "valid_png"},
    )
    assert restarted.returncode == 0, restarted.stderr
    assert (output_dir / ".complete").read_bytes() == (
        b"HRIBO REPARATION publication v1\n"
    )
    assert not backup.exists()
